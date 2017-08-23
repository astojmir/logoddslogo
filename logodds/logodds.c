/*
* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*               National Center for Biotechnology Information
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's official duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government have not placed any restriction on its use or reproduction.
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
*  Please cite the author in any work or product based on this material.
*
* ===========================================================================
*
*
*  Python wrapper to Stephen Altschul and Yi-Kuo Yu's log-odds information content codes
*
*  Wrapper author:  Aleksandar Stojmirovic
*
*
*/


#include "Python.h"

double information1(int protein, double *c, double *p);
double bildp(int dmnumber, double *c);
double bildn(double *c, double *p, double alpha);
void weighcounts(int **ords, double **counts, double *freq, int width, int numseq);
double schneider(int alpha, double *c, double *p, int flag);
extern int num_DMs;

/* Ensure that input arguments are valid */
static Py_ssize_t
check_input_sequence(PyObject *inarg, double out[20], const char *msg)
{
    Py_ssize_t out_size;
    PyObject *arr;
    PyObject *item;
    double x;
    Py_ssize_t i;

    arr = PySequence_Fast(inarg, msg);
    if (arr == NULL) {
        return 0;
    }

    out_size = PySequence_Fast_GET_SIZE(arr);
    if (out_size != 4 && out_size != 20) {
        PyErr_SetString(PyExc_TypeError, msg);
        out_size = 0;
        goto cleanup;
    }

    for (i=0; i < out_size; i++) {
        item = PySequence_Fast_GET_ITEM(arr, i);
        x = PyFloat_AsDouble(item);
        if (PyErr_Occurred()) {
            out_size = 0;
            goto cleanup;
        }
        out[i] = x;
    }
    /* Remember that arr is a new reference so it must be decremented to avoid memory leak */
 cleanup:
    Py_DECREF(arr);
    return out_size;
}

static Py_ssize_t
check_ords_row(PyObject *row, const char *msg)
{
    PyObject *arr;
    PyObject *item;
    Py_ssize_t j;
    Py_ssize_t row_size;
    long y;

    /* This checks if the row is a sequence */
    arr = PySequence_Fast(row, msg);
    if (arr == NULL) {
        row_size = 0;
        goto cleanup;
    }
    row_size = PySequence_Fast_GET_SIZE(arr);
    if (row_size < 1) {
        row_size = 0;
        PyErr_SetString(PyExc_TypeError, msg);
        goto cleanup;
    }

    /* This verifies that each entry is a long */
    /* Python integers correspond to C long type - not int. So, we are checking
       if entries are valid Python integers. Ideally, there should be more
       checks here to ensure that all entries are in the valid range for int,
       or even in a more restricted range */
    for (j=0; j < row_size; j++) {
        item = PySequence_Fast_GET_ITEM(arr, j);
        y = PyInt_AsLong(item);
        if (PyErr_Occurred()) {
            row_size = 0;
            goto cleanup;
        }
    }
 cleanup:
    Py_DECREF(arr);
    return row_size;
}


static Py_ssize_t
check_ords(PyObject *inarg, Py_ssize_t *width, const char *msg)
{
    /* Check if inarg can be converted into ords matrix. Return the number of sequences */

    /* inarg (as Python object)  must be a sequence of sequences (list of
     lists), where
     - for each i, inarg[i] have same sizes, and
     - each inarg[i][j] can be converted to int
    */

    Py_ssize_t numseq;
    Py_ssize_t current_width;
    PyObject *arr;
    PyObject *row;
    Py_ssize_t i;

    arr = PySequence_Fast(inarg, msg);
    if (arr == NULL) {
        return 0;
    }
    numseq = PySequence_Fast_GET_SIZE(arr);
    if (numseq < 1) {
        PyErr_SetString(PyExc_TypeError, msg);
        /* This is unnecessary but I left it here so we can be sure that the
           value of numseq here is exactly 0 */
        numseq = 0;
        goto cleanup;
    }

    /* First pass: check that all conditions are satisfied */
    /* I made sure that an exception is set whenever encountering an error
       during the check, and that all temporary objects are cleaned up. */
    row = PySequence_Fast_GET_ITEM(arr, 0);
    *width = check_ords_row(row, msg);
    if (*width == 0) {
        numseq = 0;
        goto cleanup;
    }
    for (i=1; i < numseq; i++) {
        row = PySequence_Fast_GET_ITEM(arr, i);
        current_width = check_ords_row(row, msg);
        if (current_width != *width) {
            /* Remember that current_width == 0 means that an error occurred
               and that the exception is already set. */
            if (!PyErr_Occurred()) {
                PyErr_SetString(PyExc_TypeError, msg);
            }
            numseq = 0;
            goto cleanup;
        }
    }

 cleanup:
    Py_DECREF(arr);
    return numseq;
}


static PyObject *
wrap_info1(PyObject *self, PyObject *args)
{
    int protein;
    double c[20];
    double p[20];
    double info;

    PyObject *obj_c;
    PyObject *obj_p;
    Py_ssize_t c_size;
    Py_ssize_t p_size;
    const char *msg_c = "counts must be a sequence of floats with size exactly 4 or 20";
    const char *msg_p = "freqs must be a sequence of floats with size exactly 4 or 20";

    if (! PyArg_ParseTuple(args, "OO", &obj_c, &obj_p)) {
        return NULL;
    }

    c_size = check_input_sequence(obj_c, c, msg_c);
    if (!c_size) {
        return NULL;
    }
    p_size = check_input_sequence(obj_p, p, msg_p);
    if (!p_size) {
        return NULL;
    }
    if (c_size != p_size) {
        PyErr_SetString(PyExc_TypeError, "counts and freqs must have the same size");
        return NULL;
    }

    protein = (c_size == 4) ? 0 : 1;
    info =  information1(protein, c, p);
    return PyFloat_FromDouble(info);
}


static PyObject *
wrap_bildp(PyObject *self, PyObject *args)
{
    int dmnumber;
    double c[20];
    double info;

    PyObject *obj_c;
    Py_ssize_t c_size;
    const char *msg_c = "counts must be a sequence of floats with size exactly 20";

    if (! PyArg_ParseTuple(args, "iO", &dmnumber, &obj_c)) {
        return NULL;
    }
    if (dmnumber < 0 || dmnumber >= num_DMs) {
        PyErr_SetString(PyExc_IndexError, "invalid dmnumber passed");
        return NULL;
    }
    c_size = check_input_sequence(obj_c, c, msg_c);
    if (!c_size) {
        return NULL;
    }
    if (c_size != 20) {
        /* check_input_sequence allows size 4 but this only covers proteins */
        PyErr_SetString(PyExc_TypeError, msg_c);
        return NULL;
    }

    info = bildp(dmnumber, c);
    return PyFloat_FromDouble(info);
}


static PyObject *
wrap_bildn(PyObject *self, PyObject *args)
{
    double c[20];
    double p[20];
    double alpha;
    double info;

    PyObject *obj_c;
    PyObject *obj_p;
    Py_ssize_t c_size;
    Py_ssize_t p_size;
    const char *msg_c = "counts must be a sequence of floats with size exactly 4";
    const char *msg_p = "freqs must be a sequence of floats with size exactly 4";

    if (! PyArg_ParseTuple(args, "OOd", &obj_c, &obj_p, &alpha)) {
        return NULL;
    }

    c_size = check_input_sequence(obj_c, c, msg_c);
    if (!c_size) {
        return NULL;
    }
    p_size = check_input_sequence(obj_p, p, msg_p);
    if (!p_size) {
        return NULL;
    }
    if (c_size != p_size || c_size != 4) {
        PyErr_SetString(PyExc_TypeError, "counts and freqs must have size 4");
        return NULL;
    }
    info = bildn(c, p, alpha);
    return PyFloat_FromDouble(info);
}


static PyObject *
wrap_weighcounts(PyObject *self, PyObject *args)
{
    /* Assume that both ords (input) and counts (output) are lists of lists */

    int **ords = NULL;
    double **counts = NULL;
    double p[20];
    Py_ssize_t width = 0;
    Py_ssize_t numseq = 0;

    PyObject *obj_ords;
    PyObject *obj_p;
    Py_ssize_t p_size;
    PyObject *obj_counts = NULL;

    PyObject *arr_ords;
    PyObject *row;
    PyObject *arr_row;
    PyObject *item;
    long y;
    Py_ssize_t i, j, k;

    const char *msg_ords = "ords must be a sequence of sequences of equal sizes";
    const char *msg_p = "freqs must be a sequence of floats with size exactly 20";

    if (! PyArg_ParseTuple(args, "OO", &obj_ords, &obj_p)) {
        return NULL;
    }

    /* Check and parse freqs */
    p_size = check_input_sequence(obj_p, p, msg_p);
    if (!p_size) {
        return NULL;
    }
    if (p_size != 20) {
        /* check_input_sequence allows size 4 but this only covers proteins */
        PyErr_SetString(PyExc_TypeError, msg_p);
        return NULL;
    }

    /* Check ords */
    numseq = check_ords(obj_ords, &width, msg_ords);
    if (!numseq) {
        return NULL;
    }

    /* Allocate ords */
    ords = (int **) calloc(sizeof(int *), numseq);
    if (ords == NULL) {
        PyErr_SetString(PyExc_MemoryError, "Could not allocate ords matrix.");
        return NULL;
    }
    for (i=0; i < numseq; i++) {
        ords[i] = (int *) malloc(width * sizeof(int));
        if (ords[i] == NULL) {
            PyErr_SetString(PyExc_MemoryError, "Could not allocate ords matrix.");
            goto cleanup_ords;
        }
    }

    /* Allocate counts */
    counts = (double **) calloc(sizeof(double *), width);
    if (counts == NULL) {
        PyErr_SetString(PyExc_MemoryError, "Could not allocate counts matrix.");
        goto cleanup_ords;
    }
    for (j=0; j < width; j++) {
        counts[j] = (double *) malloc(20 * sizeof(double));
        if (counts[j] == NULL) {
            PyErr_SetString(PyExc_MemoryError, "Could not allocate counts matrix.");
            goto cleanup_counts;
        }
    }

    /* Copy values for ords */
    /* Here we don't need to check for errors - we did that earlier */
    arr_ords = PySequence_Fast(obj_ords, msg_ords);
    for (i=0; i < numseq; i++) {
        row = PySequence_Fast_GET_ITEM(arr_ords, i);
        arr_row = PySequence_Fast(row, msg_ords);
        for (j=0; j < width; j++) {
            item = PySequence_Fast_GET_ITEM(arr_row, j);
            y = PyInt_AsLong(item);
            ords[i][j] = (int) y;
        }
        Py_DECREF(arr_row);
    }
    Py_DECREF(arr_ords);

    /* Run weighcounts (hopefully no errors can leak back) */
    weighcounts(ords, counts, p, (int) width, (int) numseq);

    /* Create a Python object to hold counts and copy values */
    obj_counts = PyList_New(width);
    if (obj_counts == NULL) {
        goto cleanup_counts;
    }
    for (j=0; j < width; j++) {
        row = PyList_New(20);
        if (row == NULL) {
            goto error_objcounts;
        }
        for (k=0; k < 20; k++) {
            item = PyFloat_FromDouble(counts[j][k]);
            if (item == NULL) {
                goto error_objcounts;
            }
            PyList_SET_ITEM(row, k, item);
        }
        PyList_SET_ITEM(obj_counts, j, row);
    }
    goto cleanup_counts;


    /* Cleanup */
 error_objcounts:
    for (j=0; j < width; j++) {
        row = PyList_GET_ITEM(obj_counts, j);
        if (row != NULL) {
            PyList_SET_ITEM(obj_counts, j, NULL);
            Py_DECREF(row);
        }
    }
    Py_DECREF(obj_counts);
    obj_counts = NULL;

 cleanup_counts:
    for (j=0; j < width; j++) {
        free(counts[j]);
    }
    free(counts);

 cleanup_ords:
    for (i=0; i < numseq; i++) {
        free(ords[i]);
    }
    free(ords);

    return obj_counts;
}


static PyObject *
wrap_schneider(PyObject *self, PyObject *args)
{
    double c[20];
    double p[20];
    int alpha;
    int flag;
    double info;

    PyObject *obj_c;
    PyObject *obj_p;
    Py_ssize_t c_size;
    Py_ssize_t p_size;
    const char *msg_c = "counts must be a sequence of floats with size exactly 4 or 20";
    const char *msg_p = "freqs must be a sequence of floats with size exactly 4 or 20";

    if (! PyArg_ParseTuple(args, "OOi", &obj_c, &obj_p, &flag)) {
        return NULL;
    }

    c_size = check_input_sequence(obj_c, c, msg_c);
    if (!c_size) {
        return NULL;
    }
    p_size = check_input_sequence(obj_p, p, msg_p);
    if (!p_size) {
        return NULL;
    }
    if (c_size != p_size ) {
        PyErr_SetString(PyExc_TypeError, "counts and freqs must have the same size");
        return NULL;
    }

    /* In general, this may be a problem, but here c_size is either 4 or 20 */
    alpha = (int) c_size;
    info = schneider(alpha, c, p, flag);

    return PyFloat_FromDouble(info);
}



PyDoc_STRVAR(info1_doc,
"_logodds.info1(counts, freqs) - Maximum-likelihood log-odds information content");
PyDoc_STRVAR(bildp_doc,
"_logodds.bildp(dmnumber, counts) - Protein-BILD information content");
PyDoc_STRVAR(bildn_doc,
"_logodds.bildn(counts, freqs, alpha) - DNA-BILD information content");
PyDoc_STRVAR(weighcounts_doc,
"_logodds.weighcounts(aligns, freqs) - weigh by effective number of observations");
PyDoc_STRVAR(schneider_doc,
"_logodds.schneider(counts, freqs, flag) - Schneider's information content");


/* Module initialization */
static PyMethodDef PyLogOdds_methods[] = {
    {"info1",  wrap_info1, METH_VARARGS, info1_doc},
    {"bildp",  wrap_bildp, METH_VARARGS, bildp_doc},
    {"bildn",  wrap_bildn, METH_VARARGS, bildn_doc},
    {"weighcounts",  wrap_weighcounts, METH_VARARGS, weighcounts_doc},
    {"schneider",  wrap_schneider, METH_VARARGS, schneider_doc},
    {NULL}  /* Sentinel */
};

PyDoc_STRVAR(PyLogOdds_doc,
"log-odds information content functions.");


#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif
PyMODINIT_FUNC init_logodds(void)
{
    PyObject* module;

    module = Py_InitModule3("_logodds", PyLogOdds_methods, PyLogOdds_doc);

    if (!module) {
        return;
    }

    PyModule_AddIntConstant(module, "num_dmixtures", num_DMs);
}

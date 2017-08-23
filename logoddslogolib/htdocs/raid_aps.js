(function(){
  window['rdv'] = {};
  var flags = 0;
  var masks = 0; 
// Validator  
  function rdvsubmit() {
   var e, m, s, f;
    if( !!(f = document.getElementById( 'irdv' )) ) {
      f.method = "GET";
      m = s = "";
      if( !!(e = document.getElementById( "isfkf" )) ) {
        if( e.value == "" ) m += "No spectrum file was uploaded.\n";
        else if( !e.value.match( /^[O,G,N,H,X,D,P]/i ) ) m += "Bad format of uploaded spectrum file.\n";
        s = e.value;
      }
      if( !!(e = document.getElementById( "isqkf" )) ) {
        if( e.value == "" ) m += "No file with peptides was uploaded.\n";
        else if( !e.value.match( /^[U,S,M,O,G,N,H,o,g,n,h]/ ) ) m += "Bad format of uploaded sequence file\n";
      }
      if( !!(e = document.getElementById( "hist" )) && e.value == "2" && s.match( /^[a-z]/ ) ) 
        m += "RAId_aPS 'Generate Histogram' mode deals only with single spectrum.\n";
      if( !!(e = document.getElementById( "db" )) && !e.value.match( /\S/ ) ) m += "No database was selected.\n";
      if( m != "" ) {
        alert( "ERROR: " + m + "\nPlease, fix the error and press 'Submit' button again." );
        return false;
      }
      if( !!(s = ao.getinput ( f, "irdvinput"  )) ) ao.setCookie( "irdvinput",  s, 366 ); else ao.delCookie( "irdvinput"  );
      if( !!(s = ao.getselect( f, "irdvselect" )) ) ao.setCookie( "irdvselect", s, 366 ); else ao.delCookie( "irdvselect" );
      if( !!(s = document.getElementById( 'irdvs' )) ) s.disabled = true;
      ao.delCookie( "tnm" );
    }
//  f.target = "_blank";
    f.submit();
    return false;
  };
//Init
  window['rdv']['initcookie'] = function() {
   var f;
    if( !!(f = ao._( 'irdv' )) ) {
      ao.setinput ( f, ao.getCookie( "irdvinput"  ) );
      ao.setselect( f, ao.getCookie( "irdvselect" ) );
    }
    ao.setEvent( 'irdv', 'submit', rdvsubmit ); 
//    ao.setEvent( 'irdv', 'reset',  rdvreset );
    if( !!(f = ao._( 'bex' )) ) {
      ao.setEvent( f,  'click',  rdvexample );
      f.disabled = false;
    }
    if( !!(f = ao._( 'irdvr' )) ) ao.setEvent( f,  'click',  rdvretrieve );
};
//Init Select
  window['rdv']['initselect'] = function() { 
   var a, i, r;
    r = RegExp( "Toggle" );
    a = document.getElementsByTagName( 'select' );
    if( a ) for( i=0; i<a.length; i++ ) {
      if( a[i].id.match( r ) ) {
        ao.updateselect( a[i] );
        ao.addEvent( a[i], 'change', function( e ) { ao.toggleselect( e ); } );
      }
    }
  };
// Update Toggle summary
  function rdvupdatesum( e ) {
   var a, i, k, f;
    f = 0;
    a = e.getElementsByTagName( 'select' );
    for( i=0; i<a.length; i++ ) {
      for( k=0; k<a[i].options.length; k++ ) {
        if( a[i].options[k].selected != a[i].options[k].defaultSelected ) f++;
      }
    }
    a = e.getElementsByTagName( 'input' );
    for( i=0; i<a.length; i++ ) {
      if( a[i].type == 'text' ) { 
        if( a[i].value != a[i].defaultValue     ) f++;
      } else if( a[i].type == 'checkbox' ) {
        if( a[i].checked != a[i].defaultChecked ) f++;
      }
    }
    return f;
  }
//Init Toggle
  window['rdv']['updatetogglesum'] = function( r ) { 
   var e, f;
    if( !!(f = ao._( r )) ) f = rdvupdatesum( f );
    if( !!(e = ao._( r + "Sum" )) ) {
      f = f ? "User settings in effect." : "Default settings used.";
      f = f + " <a href=\"#\" class=\"cas\" onclick=\"ao.toggleon(\'" + r + "Toggle\')\">Click to change</a>";
      e.innerHTML = f;
    }
  }
  window['rdv']['inittoggle'] = function() { 
   var a, i, f;
   var r = RegExp( "(Help|Toggle)$" );
    a = document.getElementsByTagName( 'a' );
    if( a ) for( i=0; i<a.length; i++ ) {
      if( a[i].id.match( r ) ) ao.addEvent( a[i], 'click', function ( e ) { //a
          e = e || window.event;
          if( e.preventDefault ) e.preventDefault();
          else if( window.event ) window.event.returnValue = false;
          ao.toggle( e );
          if( (e || window.event) == null ) return;
          e = e.srcElement ? e.srcElement : e.target;
          while( e && e.nodeType != 1 ) e = e.parentNode;
          e = e.id;
          e = e.replace( RegExp( "Toggle$" ), "" ); 
          rdv.updatetogglesum( e );
      });
    }
  };
//Reset button
  window['rdv']['reset'] = function( n ) { 
   var a, i, r, f, o;
    if(  n ) n = ao._( o = n );
    if( !n ) n = document;
    a = n.getElementsByTagName( 'select' );
    if( a ) for( i=0; i<a.length; i++ ) {
        for( r=0; r<a[i].options.length; r++ ) a[i].options[r].selected = a[i].options[r].defaultSelected;
        ao.updateselect( a[i] );
    }
    a = n.getElementsByTagName( 'input' );
    if( a ) for( i=0; i<a.length; i++ ) {
      if( a[i].type == 'checkbox' ) a[i].checked = a[i].defaultChecked;
      else if( a[i].defaultValue  ) a[i].value   = a[i].defaultValue;
    }
    if( !!(f = document.getElementById( 'irdv' )) ) {
      if( !!(a = ao.definput ( f, "irdvinput"  )) ) ao.setCookie( "irdvinput",  a, 366 ); else ao.delCookie( "irdvinput"  );
      if( !!(a = ao.defselect( f, "irdvselect" )) ) ao.setCookie( "irdvselect", a, 366 ); else ao.delCookie( "irdvselect" );
    }
    if( !!(a = ao._( "d_aadlg"  )) ) aa_saved = ao.getinput( a );
    if( !!(a = ao._( "aasmry"   )) ) a.innerHTML = updatesummary( ao._( "d_aadlg" ), "amino acids" );
    if( !!(a = ao._( "d_ptmdlg" )) ) ptm_saved = ao.getinput( a );
    if( !!(a = ao._( "ptmsmry"  )) ) a.innerHTML = updatesummary( ao._( "d_ptmdlg" ), "PTMs" );
    rdv.updatetogglesum( o );
    return true;
  };
  function rdvretrieve() {
   var e;
     if( !!(e = document.getElementById( "irdv" )) && !!e.action ) {
       e = e.action + "?ex=5";
       if( window.location ) window.location = e;
       else if( document.location.href ) document.location.href = e;
       else alert( "The option is not yet implemented: 1" );
     } else {
       alert( "The option is not yet implemented: 2" );
     }
  };
  function rdvexample() {
   var a, e, i, r, f, k;
    rdv.reset( 'mopa' );
    rdv.reset( 'sapt' );
    if( !!(a = ao._( "dsvToggle" )) ) { a.options[2].selected = true; ao.updateselect( a ); }
    if( !!(a = ao._( "d_aadlg" )) ) aa_saved = ao.getinput( a );
    if( !!(a = ao._( "aasmry"  )) ) a.innerHTML = updatesummary( ao._( "d_aadlg" ) , "amino acids" );
    if( !!(a = ao._( "d_ptmdlg")) ) ptm_saved = ao.getinput( a );
    if( !!(a = ao._( "ptmsmry" )) ) a.innerHTML = updatesummary( ao._( "d_ptmdlg" ), "PTMs" );
    k = "Dexample";
    f = "Example.dta";
    if( !!(e = document.getElementById( "isfsmry" )) ) {
      e.innerHTML = "Example spectrum uploaded.";
      e.title = "Uploaded: " + f;
    }
    if( !!(e = document.getElementById( "db" )) ) {
      e.value = "hsa";
    }
    flags = flags | 3;
    if( !!(e = document.getElementById( 'irdvs'  )) ) { e.disabled = false; e.method = "GET"; }
    if( !!(e = document.getElementById( "isfKey" )) ) e.value = "Upload";
    if( !!(e = document.getElementById( "isfkf"  )) ) e.value = k;
    if( !!(e = document.getElementById( "isfbtn" )) ) { e.innerHTML = " Upload "; e.style.display = "block"; }
    if( !!(e = document.getElementById( "isfcf"  )) ) e.innerHTML = f;
    if( !!(e = document.getElementById( "isfclr" )) ) e.style.display = "inline";
    k = "Uexample";
    f = "Example.seq";
    if( !!(e = document.getElementById( "isqsmry" )) ) {
        e.innerHTML = "Example peptide list uploaded.";
        e.title = "Uploaded: " + f;
    }
    if( !!(e = document.getElementById( "isqKey" )) ) e.value = "Upload";
    if( !!(e = document.getElementById( "isqkf"  )) ) e.value = k;
    if( !!(e = document.getElementById( "isqbtn" )) ) { e.innerHTML = " Upload "; e.style.display = "block"; }
    if( !!(e = document.getElementById( "isqcf"  )) ) e.innerHTML = f;
    if( !!(e = document.getElementById( "isqclr" )) ) e.style.display = "inline";
    rdv.updatetogglesum( 'mopa' );
    rdv.updatetogglesum( 'sapt' );
    ao.toggleon ( 'mopaToggle' );
    ao.toggleoff( 'saptToggle' );
    return true;
  };
  //Dialog
var h_init = 0;
var aa_saved = "";
var ptm_saved = "";
    function startdialog( name ) {
     var e, c;
      e = document.getElementById( "d_" + name + "dlg" );
      if( !e ) return true;
      if( h_init == 0 && name == "ptm" ) {
        h_init++;
        c = document.body.getBoundingClientRect();
        e.style.height = (c.bottom - c.top - 100).toString(10) + 'px';
      }
      e.style.display = "block";
      if( !!(c = document.getElementById( name + "tbl" )) ) c.style.height = (e.offsetHeight - 85).toString(10) + 'px';
      return false;
    };
    function filedelete( name, what ) {
     var e, k, req;
      if( !!(e = document.getElementById( 'irdvs'       )) ) e.disabled = true;
      if( !!(e = document.getElementById( name + "text" )) ) e.value="";
      if( !!(e = document.getElementById( name + "Key"  )) ) e.value="";
      if( !!(e = document.getElementById( name + "kf"   )) ) e.value="";
      if( !!(e = document.getElementById( name + "btn"  )) ) e.innerHTML=" Upload ";
      if( !!(e = document.getElementById( name + "smry" )) ){e.innerHTML="Nothing uploaded..."; e.title=""; }
      if( !!(e = document.getElementById( name + "cf"   )) ) e.innerHTML=what;
      if( !!(e = document.getElementById( name + "clr"  )) ) e.style.display = "none";
      k = ao.getCookie( name + "key" );
      ao.delCookie( name + "key" );
      ao.delCookie( name + "file" );
      if( !!k ) return false;
      if(      window.XMLHttpRequest ) req = new XMLHttpRequest();
      else if( window.ActiveXObject  ) req = new ActiveXObject("Microsoft.XMLHTTP");
      else                             return false;
      req.open( "GET", "../cgi-bin/Efile.cgi?ex=2&Key=" + k );
      req.send( null );
      return false;
    };
    function updatesummary( ed, what ) {
      var a, c, i, s;
       if( !ed || !what ) return "Error in UpdateSammary";
       a = ed.getElementsByTagName( "INPUT" );
       for( s="Selected",i=c=0; i<a.length; i++ ) if( a[i].checked ) s += (c++ == 0 ? ": " : ", ") + a[i].name;
       if( c < 1 || c > 10 ) {
         if(      c == 0 ) s = "No";
         else if( c == i ) s = "All";
         else if( c >  10 ) s = c;
         s = " " + s + " " + what + " are selected";
       }
       return s + ".";
    }
    function wrapstat( s ) { return "<span id=\"stat\">" + s + "</span>"; }
    function rdvclick( e ) {
     var a, s;
      e = e || window.event;
      if( e == null ) return false;
      if( e.preventDefault ) e.preventDefault();
      else if( window.event ) window.event.returnValue = false;
      for( e=e.srcElement?e.srcElement:e.target; e&&e.nodeType!=1; ) e = e.parentNode;
      if( e.id == 'aax' || e.id == 'aaClose' ) { 
        if( !!(a = ao._( "d_aadlg"  )) ) { 
          a.style.display = "none";
          ao.setinput( a, aa_saved ); 
        }
      } else if( e.id == 'aaAccept'  ) { 
        if( !!(a = ao._( "d_aadlg"  )) ) aa_saved = ao.getinput( a );
        if( !!(a = ao._( "aasmry"   )) ) a.innerHTML = updatesummary( ao._( "d_aadlg" ), "amino acids" );
      } else if( e.id == 'ptmx' || e.id == 'ptmClose' ) { 
        if( !!(a = ao._( "d_ptmdlg" )) ) {
          a.style.display = "none";
          ao.setinput( a, ptm_saved );
        }
      } else if( e.id == 'ptmAccept' ) {
        if( !!(a = ao._( "ptmsmry"  )) ) a.innerHTML = updatesummary( ao._( "d_ptmdlg" ), "PTMs" );
        if( !!(a = ao._( "d_ptmdlg" )) ) aa_saved = ao.getinput( a );
      } else if( e.id == 'isfx' || e.id == 'isfClose' ) { 
        if( !!(a = ao._( "d_isfdlg" )) ) a.style.display = "none";
      } else if( e.id == 'isfUpload' ) {
        if( !!(e = document.getElementById( 'irdvs'   )) ) e.disabled = true;
        if( !!(e = document.getElementById( "isfbtn"  )) ) e.className = "progress"; //e.style.display = "none";
        if( !!(e = document.getElementById( "isfsmry" )) ) e.innerHTML = wrapstat( "STATUS:" ) + " Uploading....";
        if( !!(a = ao._( "d_isfdlg" )) ) a.style.display = "none";
        if( !!(e = document.getElementById( "isfform" )) ) e.submit();
      } else if( e.id == 'isfClear' ) {
        filedelete( "isf", "no spectrum" );
        if( !!(a = ao._("d_isfbody")) ) a.innerHTML = a.innerHTML;
      } else if( e.id == 'isqx' || e.id == 'isqClose' ) { 
        if( !!(a = ao._( "d_isqdlg" )) ) a.style.display = "none";
      } else if( e.id == 'isqUpload' ) {
        if( !!(e = document.getElementById( 'irdvs'   )) ) e.disabled = true;
        if( !!(e = document.getElementById( "isqbtn"  )) ) e.className = "progress"; //e.style.display = "none";
        if( !!(e = document.getElementById( "isqsmry" )) ) e.innerHTML = wrapstat( "STATUS:" ) + " Uploading....";
        if( !!(a = ao._( "d_isqdlg" )) ) a.style.display = "none";
        if( !!(e = document.getElementById( "isqform" )) ) e.submit();
      } else if( e.id == 'isqClear' ) {
        filedelete( "isq", "no sequence" );
        if( !!(a = ao._("d_isqbody")) ) a.innerHTML = a.innerHTML;
      }
      return true;
    };
    window['rdv']['initdialog'] = function( name, what ) { 
     var et, eb, ed;
     var c, i, s;
     var a = new Array();
       if( !(ed = document.getElementById( "d_" + name + "dlg" )) ) return;
       if( !(et = document.getElementById( name + "title"      )) ) return;
       if( !(eb = document.getElementById( name + "body"       )) ) return;
       a = ed.getElementsByTagName( "*" );
       for( i=0; i<a.length; i++ ) {
         if( !!(c = a[i].className.match( /d_(\d+)/ )) ) ao.setEvent( a[i], "mousedown", ao.grab   );
         else if( a[i].className.match( /button/  )  ) ao.setEvent( a[i], "click",     rdvclick ); 
       }
       s = '<div class="b_b"> </div><div class="b_m"><span id="' + name + 'btn"> Change </span></div><div class="b_e"> </div>'
         + '<div class="b_s"><span id="'+ name + 'smry">' + updatesummary( eb, what ) + '</span></div>';
       if( !!(a = document.getElementById( "d_" + name + "test" )) ) a.innerHTML = et.innerHTML;
       if( !!(a = document.getElementById( "d_" + name + "body" )) ) {
         c = eb.parentNode;
         a.appendChild( c.removeChild( eb ) );
         c.innerHTML = s;
       }
       if( !!(c = document.getElementById( name + "btn" )) ) ao.setEvent( c, 'click', function() { return startdialog( name ); } );
    };
//Upload
    function startupload( name ) {
     var e;
      if( !(e = document.getElementById( "d_" + name + "dlg" )) ) return true;
      e.style.display = "block";
      return false;
    };
    window['rdv']['loaded'] = function( frame ) {
     var d, e, f, k, n;
      if( !(f = ao._( frame )) ) return true;
      d = f.contentWindow ? f.contentWindow.document.body : f.contentDocument ? f.contentDocument.document.body : f;
      if( !d || !d.innerHTML || !f.id.match( /(\S+)frame$/i ) ) return true;
      n = RegExp.$1;
      k = d.innerHTML.match( /Key:\s*(\S+)/ ) ? RegExp.$1 : "";
      f = "";
      if( !!(e = document.getElementById( n + "file" )) ) f = e.value.match( /([^\\\/]+)$/ ) ? RegExp.$1 : e.value;
      if( f == "" ) f = "copy/pasted text";
      if( k == "" ) f = "no file";
      if( !!(e = document.getElementById( n + "smry" )) ) {
        e.innerHTML = d.innerHTML.match( /(Error:|Status:)([^\.]+\.)/i ) 
                       ? wrapstat( RegExp.$1 ) + RegExp.$2 
                       : wrapstat( "ERROR:" ) + " Upload required.";
        e.title     = "Uploaded: " + f;
      }
      d = n == "isq" ? 2 : 1;
      flags = k != "" ? flags | d : flags & ~d;
      if( !!(e = document.getElementById( 'irdvs'    )) ){e.disabled = flags != masks; e.method = "GET"; }
      if( !!(e = document.getElementById( n + "Key"  )) ) e.value = k;
      if( !!(e = document.getElementById( n + "kf"   )) ) e.value = k;
      if( !!(e = document.getElementById( n + "btn"  )) ){e.innerHTML = k != "" ? " Reload " : " Upload "; e.className = ""; }// e.style.display = "block"; }
      if( !!(e = document.getElementById( n + "cf"   )) ) e.innerHTML = f;
      if( !!(e = document.getElementById( n + "clr"  )) ) e.style.display = k == "" ? "none" : "inline";
      if( k != "" ) { ao.setCookie( n + "key", k, 3 ); ao.setCookie( n + "file", f, 3 ); }
      else          { ao.delCookie( n + "key"       ); ao.delCookie( n + "file"       ); }
      return true;
    };
    window['rdv']['submit'] = function( form ) {
      alert( form.id + ".onsubmit()\n" );
    }
    window['rdv']['initupload'] = function( name, what, mask ) { 
     var a, c, e, ed, i, k, req;
      if( !(e = document.getElementById( name + 'body' )) ) return;
      if( !!(a = document.getElementById( name + "btn" )) ) a.className = "progress";
      masks |=  mask;
      flags &= ~mask;
      req = null;
      if( window.XMLHttpRequest ) { req = new XMLHttpRequest(); } else
      if( window.ActiveXObject  ) { req = new ActiveXObject("Microsoft.XMLHTTP"); }
      if( req == null ) return;
      e.innerHTML = '<div class="b_b"> </div><div class="b_m"><span id="' + name + 'btn"> Upload </span></div><div class="b_e"> </div>'
         + '<div class="b_s"><span id="'+ name + 'smry">Nothing uploaded...</span></div>';
      if( !!(e = document.getElementById( name + 'btn'   )) ) ao.setEvent( e, 'click', function() { return startupload( name ); } );
      if( !!(e = document.getElementById( name + 'frame' )) ) ao.addEvent( e, 'load',  function() { return rdv.loaded ( name + 'frame' ); } );
      if( !!(e = document.getElementById( name + 'clr'   )) ) ao.setEvent( e, 'click', function() { return filedelete ( name, what ); } );
      if( !(ed = document.getElementById( "d_" + name + "dlg" )) ) return;
      if( !!(e = document.getElementById( 'irdvs'        )) ) e.disabled = true;
      a = ed.getElementsByTagName( "*" );
      for( i=0; i<a.length; i++ ) {
        if( !!(c = a[i].className.match( /d_(\d+)/ )) ) ao.setEvent( a[i], "mousedown", ao.grab   );
        else if(   a[i].className.match( /button/  )  ) ao.setEvent( a[i], "click",     rdvclick ); 
      }
      if( !!(k = ao.getCookie( name + "key" )) ) {
        req.open( "GET", "../cgi-bin/Efile.cgi?Key=" + k );
        req.onreadystatechange = function() {
         var a, e, k, m;
          k = m = "";
          if( req.readyState == 4 ) { 
            a = wrapstat( "ERROR:" ) + " Upload requered: " + req.status;
            if( req.status == 200 || req.status == 0 ) { 
              a = req.responseText;
              if( (m = ao.getCookie( name + "file" )) == "" ) m = "copy/pasted text";
              if( a.match( /(Error:|Status:)([^\.]+\.)/i ) ) a = wrapstat( RegExp.$1 ) + RegExp.$2;
              if( req.responseText.match( /Key:\s*(\S+)/ ) ) {
                if( !!(RegExp.$1) ) {
                  ao.setCookie( name + "key",  k = RegExp.$1, 3 );
                  if( !!(e = document.getElementById( name + "Key" )) ) e.value = k;
                  if( !!(e = document.getElementById( name + "kf"  )) ) e.value = k;
                  if( !!(e = document.getElementById( 'irdvs'      )) ) e.disabled = ((flags |= mask) != masks);
                }
              }
            } 
            if( !!(e = document.getElementById( name + "btn" )) ) e.innerHTML=  k != "" ? " Reload " : " Upload ";
            if( k == "" ) { 
              ao.delCookie( name + "key" );
              ao.delCookie( name + "file" );
              m = "no file";
            }
            if( !!(e = document.getElementById( name + "smry" )) ) {
              e.innerHTML = a;
              e.title     = "Uploaded: " + m;
            }
            if( !!(e = document.getElementById( name + "cf"  )) ) e.innerHTML = m;
            if( !!(e = document.getElementById( name + "clr" )) ) e.style.display = k == "" ? "none" : "inline";
            if( !!(m) ) ao.setCookie( name + "file", m, 3 );
          }
        };
        req.send( null );
      }
    };
    window['rdv']['subloaded'] = function() {
     var e, d, i;
      if( !(e = document.getElementById( "subframe" )) ) { alert( "RAId: Internal error 0x00000001" ); return true; }
      i = e.contentWindow ? e.contentWindow.document : e.contentDocument ? e.contentDocument.document : null;
      d = i ? i.body : e;
      if( !d || !d.innerHTML ) { alert( "RAId: Internal error 0x00000002" ); return true; }
      e = d.innerHTML;
      d = "";
      if( !!e.match( /action="([^"]+)"/i  ) || !!e.match( /action=([^\ \>]+)/i ) ) {
        d = e = RegExp.$1;
        if( d.match( /http:\/\/[^\/]+(.+)/ ) ) e = RegExp.$1;
        window.location = e;
      } else if( !!i.getElementById && !!(e = i.getElementById( "tnm" )) ) {
        window.location = "../../qmbp-bin/Efile?ex=5&Key=" + e.innerHTML;
      }
      // else make iframe a window
      return true;
    };
    ao.addOnLoad( function() { 
      rdv.initupload( "isf", "no spectrum", 1 );
      rdv.initupload( "isq", "no sequence", 2 );
      rdv.initcookie();
      rdv.initselect();
      rdv.inittoggle();
      rdv.initdialog( "ptm", "PTMs" );
      rdv.initdialog( "aa", "amino acids" );
      rdv.updatetogglesum( "mopa" );
      rdv.updatetogglesum( "sapt" );
    } );
})();

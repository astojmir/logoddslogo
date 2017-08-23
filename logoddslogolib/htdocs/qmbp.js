(function() {
 var array_on_load = new Array();
  window.qmbp = {};
// addOnLoad event
  window.qmbp.addOnLoad = function( func ) { array_on_load.push( func ); };
  window.qmbp.onload = function() { 
    for( var i=0; i<array_on_load.length; i++ ) array_on_load[i]();
  };
// addEvent 
  window.qmbp.addEvent = function( e, evType, fn ) {
    if( typeof e == 'string' && e != "" ) e = document.getElementById( e );
    if( !e ) return false;
    if( e.addEventListener ) e.addEventListener( evType, fn, false );
    if( e.attachEvent      ) return e.attachEvent( "on" + evType, function() { return fn( window.event ); } );
    else                     e['on' + evType] = fn;
    return true;
  };
// Adjust emails. On web page it has to be in the format below:
// <a class="cat">mailuser<span class="cat"> &lt;at&gt; </span>ncbi.nlm.nih.gov</a>
  window.qmbp.adjust_emails = function() {
   var a, i;
    a = document.getElementsByTagName( "span" );
    for( i=0; i<a.length; i++ ) if( a[i].className == "cat" ) a[i].innerHTML = "@";
    a = document.getElementsByTagName( "a" );
    for( i=0; i<a.length; i++ ) if( a[i].className == "cat" ) a[i].href = "mailto:" + a[i].innerHTML.replace( /<(|\/)span[^>]*>/gi, "" );
  }
// Menu handler
  window.qmbp.initgtgl = function() {
  var a, i, f;
    a = document.getElementsByTagName( 'a' );
    if( a ) for( i=0; i<a.length; i++ ) {
      if( a[i].id.match( /^gtgl/ ) ) 
        qmbp.addEvent( a[i], 'click', function ( e ) { //a
         var ea, id, dv;
          if( (e = e || window.event) == null ) return false;
          if( e.preventDefault ) e.preventDefault();
          else if( window.event ) window.event.returnValue = false;
          e = e.srcElement ? e.srcElement : e.target;
          while( e && e.nodeType != 1 ) e = e.parentNode;
          if( !e.id.match( /^gtgl(.*)$/ ) ) return false;
          id = RegExp.$1;
          id = id.replace( RegExp( "_$" ), "" );
          if( !!(ea = document.getElementsByTagName( 'div' )) ) for( i=0; i<ea.length; i++ ) {
            if( ea[i].id.match( /^gtgl(.*)(\w)$/ ) ) {
              dv = RegExp.$1;
              if( RegExp.$2 == "a" ) {
                if( id.match( RegExp( "^" + RegExp.$1 ) ) ) 
                  ea[i].className = dv == id && ea[i].className.match( RegExp( "expanded" ) )
                                  ? ea[i].className.replace( RegExp( "expanded" ), "collapsed" )
                                  : ea[i].className.replace( RegExp( "collapsed" ), "expanded" );
                else
                  ea[i].className = ea[i].className.replace( RegExp( "expanded" ), "collapsed" );
              } else if( RegExp.$2 == "b" ) {
                if( id.match( RegExp( "^" + RegExp.$1 ) ) ) 
                  ea[i].className = dv == id && ea[i].className.match( RegExp( "opened" ) )
                                  ? ea[i].className.replace( RegExp( "opened" ), "closed" )
                                  : ea[i].className.replace( RegExp( "closed" ), "opened" );
                else                                        
                  ea[i].className = ea[i].className.replace( RegExp( "opened" ), "closed" );
              }
            }
          }
          return true;
      });
    }
  };
// Comment the line below if page does not have e-mail addresses.
  if( !! window.onload ) qmbp.addOnLoad( window.onload );
  qmbp.addOnLoad( qmbp.adjust_emails );
// Comment the line below if page does not have left menu panel
  qmbp.addOnLoad( qmbp.initgtgl );
})();
window.onload = qmbp.onload;

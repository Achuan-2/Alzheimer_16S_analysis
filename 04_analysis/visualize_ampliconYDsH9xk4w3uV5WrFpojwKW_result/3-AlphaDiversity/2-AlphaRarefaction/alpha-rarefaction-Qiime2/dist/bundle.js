webpackJsonp([0],[function(t,e,r){"use strict";function n(t){return t&&t.__esModule?t:{default:t}}var a=r(1),l=n(a);(0,l.default)()},function(t,e,r){"use strict";function n(t){return t&&t.__esModule?t:{default:t}}function a(){var t=metrics[0],e=columns[0],r=(0,l.select)("#main"),n=r.select(".controls"),a=r.select(".plotSvg"),i=r.select(".legendBoxSvg"),u=r.select(".legendTitle");c.default.initialize(t,e,n,a,i,u),(0,o.addMetricPicker)(n,metrics,t),columns.length>0?(0,o.addColumnPicker)(n,columns,e):n.selectAll(".columnPicker").remove()}Object.defineProperty(e,"__esModule",{value:!0}),e.default=a;var l=r(2),i=r(3),c=n(i),o=r(8)},,function(t,e,r){"use strict";function n(t){return t&&t.__esModule?t:{default:t}}function a(t,e){if(!(t instanceof e))throw new TypeError("Cannot call a class as a function")}function l(t,e,r,n,a,l){n.attr("href",t+".csv");var i=d[t];e&&(i=d[t][e]);var o=(0,c.setupData)(i,t);(0,u.default)(r,o,e,a,l)}Object.defineProperty(e,"__esModule",{value:!0}),e.default=void 0;var i=function(){function t(t,e){for(var r=0;r<e.length;r++){var n=e[r];n.enumerable=n.enumerable||!1,n.configurable=!0,"value"in n&&(n.writable=!0),Object.defineProperty(t,n.key,n)}}return function(e,r,n){return r&&t(e.prototype,r),n&&t(e,n),e}}(),c=r(4),o=r(5),u=n(o),s=function(){function t(){a(this,t),this.column="",this.metric="",this.svg=null,this.href=null,this.legend=null,this.legendTitle=null}return i(t,[{key:"initialize",value:function(t,e,r,n,a,i){this.href=r.select(".downloadCSV a"),this.svg=n,this.metric=t,this.column=e,this.legend=a,this.legendTitle=i,l(t,e,this.svg,this.href,this.legend,this.legendTitle)}},{key:"setColumn",value:function(t){this.column=t,l(this.metric,this.column,this.svg,this.href,this.legend,this.legendTitle)}},{key:"setMetric",value:function(t){this.metric=t,l(this.metric,this.column,this.svg,this.href,this.legend,this.legendTitle)}},{key:"getColumn",value:function(){return this.column}},{key:"getMetric",value:function(){return this.metric}},{key:"getSvg",value:function(){return this.svg}}]),t}(),f=new s;e.default=f},function(t,e,r){"use strict";function n(t){return t&&t.__esModule?t:{default:t}}function a(t,e){var r="Sequencing Depth",n=e,a=1/0,l=0,i=1/0,c=0,o=1/0,u=0,s=t.columns.indexOf("depth"),f=t.columns.indexOf("9%"),d=t.columns.indexOf("91%"),h=t.columns.indexOf("count");return t.data.forEach(function(t){var e=t[s];e<a&&(a=e),e>l&&(l=e);var r=t[f],n=t[d];r<i&&(i=r),n>c&&(c=n);var m=t[h];m>u&&(u=m),m<o&&(o=m)}),{data:t,xAxisLabel:r,yAxisLabel:n,minX:a,maxX:l,minY:i,maxY:c,minSubY:o,maxSubY:u}}function l(t,e,r){u[t]=e,u[t].dotsOpacity=1,u[t].lineOpacity=1,u[t].dots=r,u[t].line=r}function i(t,e,r){null!==e&&(u[t].dots=e,u[t].dotsOpacity="white"===e?0:1,o.default.getSvg().selectAll('[class="symbol '+t+'"]').attr("opacity",u[t].dotsOpacity)),null!==r&&(u[t].line=r,u[t].lineOpacity="white"===r?0:1,o.default.getSvg().selectAll('[class="line '+t+'"]').attr("opacity",u[t].lineOpacity))}Object.defineProperty(e,"__esModule",{value:!0}),e.curData=void 0,e.setupData=a,e.appendSeries=l,e.toggle=i;var c=r(3),o=n(c),u={};e.curData=u},function(t,e,r){"use strict";function n(t){return t&&t.__esModule?t:{default:t}}function a(t,e,r,n,a,l,o,u){function d(t,e,n){t.transition().attr("class",function(t){return"symbol "+t[A]}).attr("fill",function(t){return L(t[A])}).attr("opacity",function(t){return f.curData[t[A]].dotsOpacity}).attr("stroke",function(t){return L(t[A])}).attr("cx",function(t){return r(t[x])}).attr("cy",function(t){return n(t[e])})}function h(t,e){return(0,c.line)().x(function(t){return r(t[x])}).y(function(r){return e(r[t])})}function m(t,e,r){t.attr("class",function(t){return"line "+t.key}).attr("stroke",function(t){return L(t.key)}).attr("opacity",function(t){return f.curData[t.key].lineOpacity}).attr("fill","none").attr("d",function(t){return h(r,e)(t.values)})}var v=t.select("#chart"),p=t.select("#subChart"),y=(0,c.select)(o.node().parentNode),x=e.data.columns.indexOf("depth"),g=e.data.columns.indexOf("9%"),b=e.data.columns.indexOf("25%"),k=e.data.columns.indexOf("50%"),O=e.data.columns.indexOf("75%"),w=e.data.columns.indexOf("91%"),_=e.data.columns.indexOf("count"),A=e.data.columns.indexOf("sample-id");A===-1&&(A=e.data.columns.indexOf(l));var M=2,S=10,D=(r(e.maxX)-r(e.minX))/steps;D=D>S?S:D;var C=[e.data.data][0],P=new Set(Array.from(C,function(t){return t[A]})),L=(0,c.scaleOrdinal)(c.schemeCategory20).domain(P),j=Array.from(P);o.selectAll(".legend").remove(),u.selectAll(".legend").remove(),o.attr("height",20*j.length);var Y=0,N="Select%20All";(0,f.appendSeries)(N,[],"black"),(0,f.toggle)(N,null,null),(0,s.default)(u,N,10,L);var T=j.sort(function(t,e){var r=isNaN(t),n=isNaN(e);return r&&n?t>e?1:-1:!r&&n?1:r&&!n?-1:t-e}),X=!0,I=!1,E=void 0;try{for(var z,B=function(){var t=i(z.value,2),e=t[0],r=t[1];Y=20*(e+.5);var n=C.filter(function(t){return t[A]===r}).sort(function(t,e){return t[x]-e[x]}),a=L(r);(0,f.appendSeries)(r,n,a),(0,f.toggle)(r,null,null),(0,s.default)(o,r,Y,L)},q=T.entries()[Symbol.iterator]();!(X=(z=q.next()).done);X=!0)B()}catch(t){I=!0,E=t}finally{try{!X&&q.return&&q.return()}finally{if(I)throw E}}var J=p.selectAll(".symbol").data(C);J.exit().remove();var R=J.enter().append("circle").attr("r",M);J.call(d,_,a),R.call(d,_,a),y.attr("height",""+(Y+10)).attr("width","200");var U=(0,c.nest)().key(function(t){return t[A]}).entries(C),V=v.selectAll(".line").data(U);V.exit().remove(),V.enter().append("path").call(m,n,k),V.call(m,n,k);var F=v.selectAll(".symbol").data(e.data.data);F.exit().remove();var G=F.enter().append("g"),H=F.merge(G).attr("class",function(t){return"symbol "+t[A]}).attr("opacity",function(t){return f.curData[t[A]].dotsOpacity}).attr("transform",function(t){return"translate("+r(t[x])+", 0)"}),K=H.selectAll("line.center").data(function(t){return[t]});K.exit().remove();var Q=K.enter().append("line");K.merge(Q).attr("class","center").attr("x1",0).attr("y1",function(t){return n(t[g])}).attr("x2",0).attr("y2",function(t){return n(t[w])}).attr("stroke-width",1).attr("stroke",function(t){return L(t[A])});var W=H.selectAll("rect.box").data(function(t){return[t]});W.exit().remove();var Z=W.enter().append("rect");W.merge(Z).attr("class","box").attr("x",-(D/2)).attr("y",function(t){return n(t[O])}).attr("width",D).attr("height",function(t){return n(t[b])-n(t[O])}).attr("fill","white").attr("stroke-width",1).attr("stroke",function(t){return L(t[A])});var $=H.selectAll("line.median").data(function(t){return[t]});$.exit().remove();var tt=$.enter().append("line");$.merge(tt).attr("class","median").attr("x1",-(D/2)).attr("y1",function(t){return n(t[k])}).attr("x2",D/2).attr("y2",function(t){return n(t[k])}).attr("stroke-width",1).attr("stroke",function(t){return L(t[A])});var et=p.selectAll(".line").data(U);et.exit().remove(),et.enter().append("path").call(m,a,_),et.call(m,a,_)}function l(t,e,r,n,l){function i(t,e,r){var n=.03*(e-t);if(Number.isInteger(t)&&Number.isInteger(e)){n=Math.max(Math.round(n),1);var a=Math.max(3,e-t+2*n);r.ticks(Math.min(a,12),"d")}return n}var u=400,s=1e3,f={top:20,left:80,right:50,bottom:50},d=t.select("#chart"),h=t.select("#subChart"),m=e.xAxisLabel,v=e.yAxisLabel,p=e.minX,y=e.maxX,x=e.minY,g=e.maxY,b=e.minSubY,k=e.maxSubY,O=(0,c.axisBottom)(),w=(0,c.axisLeft)(),_=(0,c.axisLeft)(),A=i(p,y,O),M=i(b,k,_);M=1===M?2:M;var S=(0,c.scaleLinear)().domain([p-A,y+A]).range([0,s]).nice(),D=(0,c.scaleLinear)().domain([x,g]).range([u,0]).nice(),C=(0,c.scaleLinear)().domain([b-M,k+M]).range([u,0]).nice();O.scale(S),w.scale(D),_.scale(C),(0,o.setupXLabel)(t,s,u,m,O);var P=(0,o.setupYLabels)(t,u,v,w,_),L=Math.max(f.left,P);t.attr("width",s+L+f.right).attr("height",2*(u+f.bottom+f.top)),(0,c.select)(t.node().parentNode).style("width",s+L+f.right+"px").style("height",2*(u+f.bottom+f.top)+"px"),d.attr("transform","translate("+L+","+f.top+")"),h.attr("transform","translate("+L+","+(u+f.bottom+f.top)+")"),a(t,e,S,D,C,r,n,l)}Object.defineProperty(e,"__esModule",{value:!0});var i=function(){function t(t,e){var r=[],n=!0,a=!1,l=void 0;try{for(var i,c=t[Symbol.iterator]();!(n=(i=c.next()).done)&&(r.push(i.value),!e||r.length!==e);n=!0);}catch(t){a=!0,l=t}finally{try{!n&&c.return&&c.return()}finally{if(a)throw l}}return r}return function(e,r){if(Array.isArray(e))return e;if(Symbol.iterator in Object(e))return t(e,r);throw new TypeError("Invalid attempt to destructure non-iterable instance")}}();e.default=l;var c=r(2),o=r(6),u=r(7),s=n(u),f=r(4)},function(t,e,r){"use strict";function n(t,e,r,n,a){t.select("#chart .x.axis").attr("transform","translate(0,"+r+")").transition().call(a),t.select("#chart .x.label").attr("text-anchor","middle").style("font","12px sans-serif").text(n).attr("transform","translate("+e/2+","+(r+30)+")"),t.select("#subChart .x.axis").attr("transform","translate(0,"+r+")").transition().call(a),t.select("#subChart .x.label").attr("text-anchor","middle").style("font","12px sans-serif").text(n).attr("transform","translate("+e/2+","+(r+30)+")")}function a(t,e,r,n,a){function i(t,e){var r=t.attr("text-anchor","middle").style("font","12px sans-serif").text(e);return r}var c=t.select("#chart .y.axis").call(n),o=t.select("#chart .y.label").call(i,r),u=Array.from(c.selectAll("text")._groups[0]).map(function(t){return t.getComputedTextLength()}),s=(0,l.max)(u)+20;o.attr("transform","translate(-"+s+","+e/2+")rotate(-90)"),t.select("#subChart .y.axis").call(a);var f=t.select("#subChart .y.label").call(i,"Number of samples");return f.attr("transform","translate(-"+s+","+e/2+")rotate(-90)"),s}Object.defineProperty(e,"__esModule",{value:!0}),e.setupXLabel=n,e.setupYLabels=a;var l=r(2)},function(t,e,r){"use strict";function n(t,e,r,n){var i="Select%20All",c=t.append("rect").attr("id","rect"+e).attr("class","legend rect").attr("x",2).attr("y",r-2.5).attr("width",15).attr("height",5).attr("stroke","black").attr("stroke-width","1").attr("fill",l.curData[e].line);c.on("click",function(){var t=0===l.curData[e].lineOpacity;if(e===i){var r=!0,o=!1,u=void 0;try{for(var s,f=Object.keys(l.curData)[Symbol.iterator]();!(r=(s=f.next()).done);r=!0){var d=s.value,h=(0,a.select)('[id="rect'+d+'"]'),m=t?n(d):"white";(0,l.toggle)(d,null,m),h.attr("fill",l.curData[d].line)}}catch(t){o=!0,u=t}finally{try{!r&&f.return&&f.return()}finally{if(o)throw u}}}var v=e===i?"black":n(e),p=t?v:"white";(0,l.toggle)(e,null,p),c.attr("fill",l.curData[e].line)});var o=t.append("circle").attr("id","dot"+e).attr("class","legend circle").attr("cx",30).attr("cy",r).attr("r",5).attr("stroke","black").attr("stroke-width","1").attr("fill",l.curData[e].dots);o.on("click",function(){var t=0===l.curData[e].dotsOpacity;if(e===i){var r=!0,c=!1,u=void 0;try{for(var s,f=Object.keys(l.curData)[Symbol.iterator]();!(r=(s=f.next()).done);r=!0){var d=s.value,h=(0,a.select)('[id="dot'+d+'"]'),m=t?n(d):"white";(0,l.toggle)(d,m,null),h.attr("fill",l.curData[d].dots)}}catch(t){c=!0,u=t}finally{try{!r&&f.return&&f.return()}finally{if(c)throw u}}}var v=e===i?"black":n(e),p=t?v:"white";(0,l.toggle)(e,p,null),o.attr("fill",l.curData[e].dots)}),t.append("text").attr("class","legend").attr("x",40).attr("y",r+5).style("font","10px sans-serif").text(decodeURI(e))}Object.defineProperty(e,"__esModule",{value:!0}),e.default=n;var a=r(2),l=r(4)},function(t,e,r){"use strict";function n(t){return t&&t.__esModule?t:{default:t}}function a(t,e,r){t.select(".metricPicker select").on("change",function(){var t=e[this.selectedIndex];c.default.setMetric(t)}).selectAll("option").data(e).enter().append("option").attr("value",function(t){return t}).text(function(t){return t}).property("selected",function(t){return t===r})}function l(t,e,r){t.select(".columnPicker select").on("change",function(){var t=e[this.selectedIndex];c.default.setColumn(t)}).selectAll("option").data(e).enter().append("option").attr("value",function(t){return t}).text(function(t){return t}).property("selected",function(t){return t===r})}Object.defineProperty(e,"__esModule",{value:!0}),e.addMetricPicker=a,e.addColumnPicker=l;var i=r(3),c=n(i)}]);
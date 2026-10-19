/* Dark defaults for every Plotly chart on the site */
(function(){if(!window.Plotly)return;var o=Plotly.newPlot,ax=function(a){return Object.assign({gridcolor:"#16233d",zerolinecolor:"#3b4f6e",linecolor:"#2a3f63"},a||{})};
Plotly.newPlot=function(id,d,l,c){l=l||{};return o.call(Plotly,id,d,Object.assign({},l,{font:Object.assign({color:"#c9d6e8",family:"Latin Modern, CM Greek, CM Symbols, serif"},l.font||{}),xaxis:ax(l.xaxis),yaxis:ax(l.yaxis),paper_bgcolor:"rgba(0,0,0,0)",plot_bgcolor:"rgba(0,0,0,0)"}),c)}})();

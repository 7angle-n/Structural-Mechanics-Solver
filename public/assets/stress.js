/* Bending + shear (+ axial) stress distribution across a section. y is measured up from the centroid.
   Sagging M (+) puts the top in compression: sigma = N/A - M*y/I. tau = V*Q/(I*b). Consistent units, e.g. N, mm, MPa. */
const g=v=>(+v.toPrecision(5)).toString();
const DIM={rect:[["b",100],["h",200]],circ:[["d",120]],isec:[["bf",120],["tf",15],["tw",10],["hw",200]]};
function dims(){$("dims").innerHTML=DIM[$("sh").value].map(([n,d])=>`<label>${n} <input type="number" class="dm" value="${d}" step="any"></label>`).join("")}
function stress(){
 const sh=$("sh").value,v=[...document.querySelectorAll(".dm")].map(e=>+e.value),M=+$("M").value||0,V=+$("V").value||0,N=+$("N").value||0,pi=Math.PI;
 let A,I,c,tau,L=[{h:"Section properties"}];
 if(sh=="rect"){const[b,h]=v;A=b*h;I=b*h**3/12;c=h/2;tau=y=>3*V/(2*A)*(1-(y/c)**2);L.push(`A=${g(A)},\\quad I=\\frac{bh^3}{12}=${g(I)},\\quad c=\\frac h2=${g(c)}`)}
 if(sh=="circ"){const[d]=v;A=pi*d*d/4;I=pi*d**4/64;c=d/2;tau=y=>4*V/(3*A)*(1-(y/c)**2);L.push(`A=${g(A)},\\quad I=\\frac{\\pi d^4}{64}=${g(I)},\\quad c=${g(c)}`)}
 if(sh=="isec"){const[bf,tf,tw,hw]=v,H=hw+2*tf;A=bf*tf*2+tw*hw;I=(bf*H**3-(bf-tw)*hw**3)/12;c=H/2;
  tau=y=>{y=Math.abs(y);return y<=hw/2?V/(I*tw)*(bf*tf*(hw/2+tf/2)+tw/2*(hw*hw/4-y*y)):V/(I*bf)*(bf/2)*(c*c-y*y)};
  L.push(`H=h_w+2t_f=${g(H)},\\quad A=2b_ft_f+t_wh_w=${g(A)}`,`I=\\frac{b_fH^3-(b_f-t_w)h_w^3}{12}=${g(I)},\\quad c=${g(c)}`)}
 const S=I/c,sTop=N/A-M*c/I,sBot=N/A+M*c/I,tm=tau(0);
 L.push({h:"Bending stress"},`\\sigma(y)=\\frac{N}{A}-\\frac{My}{I},\\quad S=\\frac{I}{c}=${g(S)}`,`\\sigma_{top}=\\frac NA-\\frac{Mc}{I}=${g(sTop)},\\quad \\sigma_{bottom}=\\frac NA+\\frac{Mc}{I}=${g(sBot)}`,`\\text{(negative = compression; sagging }M\\text{ compresses the top)}`);
 if(N)L.push(`\\text{Neutral axis shifts to }y=\\frac{NI}{AM}=${M?g(N*I/(A*M)):"\\infty"}`);
 L.push({h:"Shear stress"},sh=="rect"?`\\tau(y)=\\frac{3V}{2A}\\left(1-\\frac{y^2}{c^2}\\right),\\quad \\tau_{max}=\\frac{3V}{2A}=${g(tm)}`:sh=="circ"?`\\tau(y)=\\frac{4V}{3A}\\left(1-\\frac{y^2}{c^2}\\right),\\quad \\tau_{max}=\\frac{4V}{3A}=${g(tm)}`:
  `\\text{Web: }\\tau=\\frac{V}{It_w}\\left[b_ft_f\\frac{h_w+t_f}{2}+\\frac{t_w}{2}\\left(\\frac{h_w^2}{4}-y^2\\right)\\right],\\quad \\tau_{max}=\\tau(0)=${g(tm)}`,`\\text{Flange }(h_w/2\\le|y|\\le c):\\ \\tau=\\frac{V}{I\\,b_f}\\cdot\\frac{b_f}{2}(c^2-y^2)=\\frac{V}{2I}\\left(c^2-y^2\\right)`);
 $("ans").innerHTML=`<div>I<b>${g(I)}</b></div><div>S = I/c<b>${g(S)}</b></div><div>σ top<b>${g(sTop)}</b></div><div>σ bottom<b>${g(sBot)}</b></div><div>τ max (at NA)<b>${g(tm)}</b></div>`;
 steps($("st"),L);
 const ys=Array.from({length:101},(_,i)=>c-2*c*i/100),lay=(t,x)=>({title:t,xaxis:{title:x,zeroline:true},yaxis:{title:"y from centroid"},margin:{t:36,l:60,r:12,b:40},height:330,paper_bgcolor:"rgba(0,0,0,0)",plot_bgcolor:"rgba(0,0,0,0)"});
 Plotly.newPlot("p1",[{x:ys.map(y=>N/A-M*y/I),y:ys,fill:"tozerox",line:{color:"#c0392b"}}],lay("Bending (normal) stress","σ"),{responsive:true});
 Plotly.newPlot("p2",[{x:ys.map(tau),y:ys,fill:"tozerox",line:{color:"#1d4ed8"}}],lay("Shear stress","τ"),{responsive:true});
}

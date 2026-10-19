/* Column buckling: Euler, Johnson (intermediate), secant (eccentric). Consistent units, e.g. mm, MPa, N. */
const g=v=>(+v.toPrecision(5)).toString();
function col(){
 const sh=$("sh").value,v=[...document.querySelectorAll(".dm")].map(e=>+e.value),E=+$("E").value,Lc=+$("Lc").value,K=+$("K").value,sy=+$("sy").value||0,P=+$("P").value||0,e=+$("e").value||0,FS=+$("fs").value||1,pi=Math.PI;
 let A,I,c,L=[{h:"Section properties"}];
 if(sh=="rect"){const[b,h]=v;A=b*h;I=Math.min(b*h**3,h*b**3)/12;c=Math.min(b,h)/2;L.push(`A=bh=${g(A)},\\quad I_{min}=\\frac{${b<=h?"h\\,b^3":"b\\,h^3"}}{12}=${g(I)}`)}
 if(sh=="circ"){const[d]=v;A=pi*d*d/4;I=pi*d**4/64;c=d/2;L.push(`A=\\frac{\\pi d^2}{4}=${g(A)},\\quad I=\\frac{\\pi d^4}{64}=${g(I)}`)}
 if(sh=="hcirc"){const[D,d]=v;A=pi*(D*D-d*d)/4;I=pi*(D**4-d**4)/64;c=D/2;L.push(`A=${g(A)},\\quad I=\\frac{\\pi(D^4-d^4)}{64}=${g(I)}`)}
 if(sh=="gen"){[A,I,c]=v}
 const r=Math.sqrt(I/A),Le=K*Lc,lam=Le/r,Pe=pi*pi*E*I/Le**2,lc=sy?pi*Math.sqrt(2*E/sy):Infinity,euler=lam>=lc;
 const sj=sy*(1-sy*lam*lam/(4*pi*pi*E)),Pcr=euler?Pe:sj*A,scr=Pcr/A;
 L.push(`r=\\sqrt{I/A}=${g(r)}`,{h:"Effective length and slenderness"},`L_e=KL=${K}\\times${Lc}=${g(Le)},\\quad \\lambda=\\frac{L_e}{r}=${g(lam)}`);
 if(sy)L.push(`\\lambda_c=\\pi\\sqrt{\\frac{2E}{\\sigma_y}}=${g(lc)}\\ \\Rightarrow\\ \\text{${euler?"long column: Euler applies":"intermediate column: use Johnson"}}`);
 L.push({h:euler?"Euler critical load":"Johnson critical stress"});
 if(euler)L.push(`P_{cr}=\\frac{\\pi^2EI}{(KL)^2}=\\frac{\\pi^2(${g(E)})(${g(I)})}{${g(Le*Le)}}=${g(Pcr)},\\quad \\sigma_{cr}=\\frac{P_{cr}}{A}=${g(scr)}`);
 else L.push(`\\sigma_{cr}=\\sigma_y\\left[1-\\frac{\\sigma_y}{4\\pi^2E}\\left(\\frac{L_e}{r}\\right)^2\\right]=${g(sj)},\\quad P_{cr}=\\sigma_{cr}A=${g(Pcr)}`,`\\text{(Euler would give }${g(Pe)}\\text{, unconservative here)}`);
 L.push(`P_{allow}=\\frac{P_{cr}}{FS}=\\frac{${g(Pcr)}}{${FS}}=${g(Pcr/FS)}`);
 let extra="";
 if(P){L.push({h:"Check applied load P = "+P},`P/P_{allow}=${g(P/(Pcr/FS))}\\ \\Rightarrow\\ \\text{${P<=Pcr/FS?"safe":"NOT safe"}}`);extra=`<div>${P<=Pcr/FS?"Safe":"Unsafe"}<b>${g(P/(Pcr/FS))}</b>P / P_allow</div>`}
 if(P&&e){const sm=P/A*(1+e*c/(r*r)/Math.cos(Le/(2*r)*Math.sqrt(P/(E*A))));L.push({h:"Eccentric load (secant formula)"},`\\sigma_{max}=\\frac{P}{A}\\left[1+\\frac{ec}{r^2}\\sec\\left(\\frac{L_e}{2r}\\sqrt{\\frac{P}{EA}}\\right)\\right]=${g(sm)}${sy?`\\ \\Rightarrow\\ \\text{${sm<=sy?"below":"ABOVE"} yield}`:""}`);extra+=`<div>σmax (secant)<b>${g(sm)}</b></div>`}
 $("ans").innerHTML=`<div>λ<b>${g(lam)}</b></div><div>Pcr<b>${g(Pcr)}</b>${euler?"Euler":"Johnson"}</div><div>σcr<b>${g(scr)}</b></div><div>P allowable<b>${g(Pcr/FS)}</b></div>`+extra;
 steps($("st"),L);
 const xs=Array.from({length:250},(_,i)=>i+1),eu=xs.map(l=>pi*pi*E/l/l),tr=[{x:xs,y:eu.map(y=>sy?Math.min(y,sy*1.6):y),name:"Euler",line:{dash:"dot",color:"#7d8f87"}}];
 if(sy)tr.push({x:xs.filter(l=>l<=lc),y:xs.filter(l=>l<=lc).map(l=>sy*(1-sy*l*l/(4*pi*pi*E))),name:"Johnson",line:{color:"#4f7c6a"}},{x:xs.filter(l=>l>=lc),y:eu.filter((_,i)=>xs[i]>=lc),name:"Euler (valid)",line:{color:"#1d4ed8"}});
 tr.push({x:[lam],y:[scr],mode:"markers",name:"Your column",marker:{size:12,color:"#c0392b"}});
 Plotly.newPlot("pl",tr,{title:"Critical stress vs slenderness ratio",xaxis:{title:"L_e / r",range:[0,250]},yaxis:{title:"σcr",range:[0,sy?sy*1.2:Math.min(...eu)*4]},margin:{t:36,l:64,r:12,b:40},height:300,paper_bgcolor:"rgba(0,0,0,0)",plot_bgcolor:"rgba(0,0,0,0)"},{responsive:true});
}
const DIM={rect:[["b",50],["h",100]],circ:[["d",60]],hcirc:[["D",80],["d",60]],gen:[["A",5000],["I",1e6],["c",25]]};
function dims(){$("dims").innerHTML=DIM[$("sh").value].map(([n,d])=>`<label>${n} <input type="number" class="dm" value="${d}" step="any"></label>`).join("")}

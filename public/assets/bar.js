/* Stepped bar (axial) or shaft (torsion). Same maths: member "stiffness" k = EA or GJ.
   Left end fixed; right end free or fixed. Right-fixed case uses the compatibility equation. */
const g=v=>(+v.toPrecision(5)).toString();
const stepXY=(xs,v)=>({x:xs.flatMap((x,i)=>i<v.length?[x,xs[i+1]]:[]),y:v.flatMap(a=>[a,a])});
function bar(mode){
 const ax=mode=="axial",need=ax?4:5;
 const S=$("tx").value.trim().split("\n").map(l=>l.trim().split(/\s+/).map(Number)).filter(a=>a.length>=need&&a.every(v=>!isNaN(v))).map(a=>{
  if(ax){const[L,A,E,P,al=0,dT=0]=a;return{L,P,k:A*E,A,E,th:al*dT*L}}const[L,D,d,G,T]=a,J=Math.PI*(D**4-d**4)/32;return{L,P:T,k:G*J,D,d,G,J,th:0}});
 if(!S.length){$("ans").innerHTML="<p>Enter at least one segment.</p>";return}
 const n=S.length,both=$("sup").value=="both",nm=ax?"N":"T",st=ax?"\\sigma":"\\tau",df=ax?"\\delta":"\\varphi";
 const Si=S.map((_,i)=>S.slice(i).reduce((s,r)=>s+r.P,0)),fl=S.map(r=>r.L/r.k),sf=fl.reduce((a,b)=>a+b,0);
 const thT=S.reduce((s,r)=>s+r.th,0),free=Si.reduce((s,v,i)=>s+v*fl[i],0)+thT,gap=+(($("gap")||{}).value)||0,contact=both&&free>gap,X=contact?-(free-gap)/sf:0,N=Si.map(v=>v+X),RA=-(Si[0]+X);
 const dl=N.map((v,i)=>v*fl[i]+S[i].th),u=[0];dl.forEach(v=>u.push(u[u.length-1]+v));
 const sg=N.map((v,i)=>ax?v/S[i].A:v*(S[i].D/2)/S[i].J),xs=[0];S.forEach(r=>xs.push(xs[xs.length-1]+r.L));
 const mx=Math.max(...sg.map(Math.abs));
 $("ans").innerHTML=`<div>${ax?"Left reaction":"Left reaction torque"}<b>${g(RA)}</b></div>${both?`<div>Right reaction<b>${g(X)}</b></div>`:""}<div>Max |${ax?"σ":"τ"}|<b>${g(mx)}</b></div><div>${ax?"Total elongation":"Total twist"}<b>${g(u[n])}</b>${ax?"":g(u[n]*180/Math.PI)+"°"}</div>`;
 const L=[{h:"Equilibrium: internal "+(ax?"force":"torque")+" in each segment"},ax?"\\text{Cut a section; positive }N\\text{ = tension. Sum the external loads to the right:}":"\\text{Cut a section; positive }T\\text{ by the right-hand rule. Sum the external torques to the right:}"];
 S.forEach((r,i)=>L.push(`${nm}_{${i+1}}=${S.slice(i).map(q=>g(q.P)).join("+")}${both?"+X":""}`));
 if(both){L.push({h:"Compatibility (both ends fixed)"},`\\text{Total ${ax?"elongation":"twist"} with the right end released: }\\sum\\frac{S_iL_i}{k_i}${thT?"+\\sum\\alpha\\Delta T\\,L":""}=${g(free)}`);
  if(!contact)L.push(`${g(free)}\\le\\text{gap}=${g(gap)}\\ \\Rightarrow\\ \\text{the bar never touches the right support, so }X=0`);
  else L.push(`\\sum\\frac{(S_i+X)L_i}{k_i}${thT?"+\\sum\\alpha\\Delta T\\,L":""}=${g(gap)}\\ \\Rightarrow\\ X=-\\frac{${g(free)}-${g(gap)}}{${g(sf)}}=${g(X)}`,`R_{left}=-(\\sum P+X)=${g(RA)}`,...S.map((_,k)=>`${nm}_{${k+1}}=${g(Si[k])}+(${g(X)})=${g(N[k])}`))}
 L.push({h:ax?"Stress, strain, elongation":"Shear stress and twist"});
 S.forEach((r,i)=>{L.push({h:"Segment "+(i+1)});
  if(ax)L.push(`\\sigma=\\frac{N}{A}=\\frac{${g(N[i])}}{${g(r.A)}}=${g(sg[i])},\\quad \\varepsilon=\\frac{\\sigma}{E}=${g(sg[i]/r.E)},\\quad \\delta=\\frac{NL}{AE}=${g(dl[i])}`);
  else L.push(`J=\\frac{\\pi(D^4-d^4)}{32}=${g(r.J)},\\quad \\tau_{max}=\\frac{T\\,D/2}{J}=${g(sg[i])},\\quad \\gamma=\\frac{\\tau}{G}=${g(sg[i]/r.G)}`,`\\varphi=\\frac{TL}{GJ}=${g(dl[i])}\\ \\text{rad}=${g(dl[i]*180/Math.PI)}^\\circ`)});
 L.push({h:"Node displacements"},`${ax?"u":"\\theta"}=${u.map(g).join(",\\ ")}`);
 steps($("st"),L);
 const lay=(t,y)=>({title:t,margin:{t:36,l:64,r:12,b:34},xaxis:{title:"x"},yaxis:{title:y},height:230,paper_bgcolor:"rgba(0,0,0,0)",plot_bgcolor:"rgba(0,0,0,0)"});
 const a=stepXY(xs,N),b=stepXY(xs,sg);
 Plotly.newPlot("p1",[{...a,fill:"tozeroy",line:{color:"#4f7c6a"}}],lay(ax?"Axial force diagram":"Torque diagram",nm),{responsive:true});
 Plotly.newPlot("p2",[{...b,fill:"tozeroy",line:{color:"#c0392b"}}],lay(ax?"Stress diagram (uniform in each segment)":"Max shear stress diagram",ax?"σ":"τ"),{responsive:true});
 Plotly.newPlot("p3",[{x:xs,y:u,line:{color:"#1d4ed8"}}],lay(ax?"Displacement of cross-sections":"Angle of twist of cross-sections",ax?"u":"θ (rad)"),{responsive:true});
 if(!ax){const i=sg.map(Math.abs).indexOf(mx),r=S[i],R=r.D/2,ri=r.d/2,t=N[i]/r.J,gp=(a,c)=>Array.from({length:21},(_,k)=>a+(c-a)*k/20);
  const rr=[...gp(-R,-ri),null,...gp(ri,R)];Plotly.newPlot("p4",[{x:rr,y:rr.map(v=>v===null?null:t*v),fill:"tozeroy",line:{color:"#a0612b"}}],lay("Shear stress across the section (segment "+(i+1)+"): τ = T·r/J","τ"),{responsive:true})}
}

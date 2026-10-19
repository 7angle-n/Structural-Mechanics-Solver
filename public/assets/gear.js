/* Gear pair (shaft 1 -> gears B,C -> shaft 2 fixed at D) and shaft sizing from power. SI units: m, Pa, N·m. */
const gv=id=>+$(id).value,G5=v=>(+v.toPrecision(5)).toString();
function gear(){
 const pi=Math.PI,T=gv("gT"),r1=gv("gr1"),r2=gv("gr2"),sh=n=>{const L=gv("gL"+n),D=gv("gD"+n),d=gv("gd"+n),G=gv("gG"+n);return{L,D,d,G,J:pi*(D**4-d**4)/32}},a=sh(1),b=sh(2);
 const F=T/r1,T2=F*r2,t1=T*(a.D/2)/a.J,t2=T2*(b.D/2)/b.J,pAB=T*a.L/(a.G*a.J),pC=T2*b.L/(b.G*b.J),pB=pC*r2/r1,pA=pB+pAB,dg=x=>G5(x*180/pi);
 $("ans2").innerHTML=`<div>Tooth force F<b>${G5(F)}</b></div><div>Torque in shaft 2<b>${G5(T2)}</b></div><div>τmax shaft 1<b>${G5(t1)}</b></div><div>τmax shaft 2<b>${G5(t2)}</b></div><div>Twist at A<b>${G5(pA)} rad</b>${dg(pA)}°</div>`;
 steps($("st2"),[{h:"Gear mesh"},`F=\\frac{T}{r_1}=\\frac{${G5(T)}}{${r1}}=${G5(F)},\\quad T_2=F\\,r_2=T\\frac{r_2}{r_1}=${G5(T2)}`,`\\text{Gear B and C turn in opposite senses; speed ratio }\\frac{n_2}{n_1}=\\frac{r_1}{r_2}=${G5(r1/r2)}`,
  {h:"Shear stress"},`J_1=\\frac{\\pi(D^4-d^4)}{32}=${G5(a.J)},\\quad \\tau_1=\\frac{TD/2}{J_1}=${G5(t1)}`,`J_2=${G5(b.J)},\\quad \\tau_2=\\frac{T_2D/2}{J_2}=${G5(t2)}`,
  {h:"Twist (shaft 2 is fixed at D)"},`\\varphi_C=\\frac{T_2L_2}{G_2J_2}=${G5(pC)}\\ \\text{rad}`,`\\text{Gears: }r_1\\varphi_B=r_2\\varphi_C\\ \\Rightarrow\\ \\varphi_B=\\varphi_C\\frac{r_2}{r_1}=${G5(pB)}`,
  `\\varphi_{A/B}=\\frac{TL_1}{G_1J_1}=${G5(pAB)}`,`\\varphi_A=\\varphi_B+\\varphi_{A/B}=${G5(pA)}\\ \\text{rad}=${dg(pA)}^\\circ`])}
function design(){
 const pi=Math.PI,P=gv("dP")*1000,n=gv("dn"),ta=gv("dt")*1e6,k=gv("dk")||0,G=gv("dG")||0,lim=gv("dl")||0,T=P*60/(2*pi*n);
 const dS=Math.cbrt(16*T/(pi*ta)),Do=Math.cbrt(16*T/(pi*ta*(1-k**4))),L=[{h:"Torque from power"},`T=\\frac{P}{\\omega}=\\frac{60P}{2\\pi n}=\\frac{60(${G5(P)})}{2\\pi(${n})}=${G5(T)}\\ \\text{N}\\cdot\\text{m}`,{h:"Strength (shear stress limit)"},
  `\\text{Solid: }d=\\left(\\frac{16T}{\\pi\\tau_{all}}\\right)^{1/3}=${G5(dS*1000)}\\ \\text{mm}`,`\\text{Hollow, }k=d_i/D=${k}:\\ D=\\left(\\frac{16T}{\\pi\\tau_{all}(1-k^4)}\\right)^{1/3}=${G5(Do*1000)}\\ \\text{mm},\\ d_i=${G5(k*Do*1000)}\\ \\text{mm}`];
 let ds=dS;if(G&&lim){const J=T/(G*1e9*lim*pi/180),dr=Math.pow(32*J/pi,.25);L.push({h:"Stiffness (twist per metre limit)"},`J\\ge\\frac{T}{G\\,\\theta'}=${G5(J)},\\quad d\\ge\\left(\\frac{32J}{\\pi}\\right)^{1/4}=${G5(dr*1000)}\\ \\text{mm}`);ds=Math.max(dS,dr);L.push(`d_{solid}=\\max(\\text{strength},\\text{stiffness})=${G5(ds*1000)}\\ \\text{mm}`)}
 $("ans3").innerHTML=`<div>Torque<b>${G5(T)}</b>N·m</div><div>Solid shaft d<b>${G5(ds*1000)}</b>mm</div><div>Hollow D / d<b>${G5(Do*1000)} / ${G5(k*Do*1000)}</b>mm</div>`;steps($("st3"),L)}

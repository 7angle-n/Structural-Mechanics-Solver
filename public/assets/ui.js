const $=id=>document.getElementById(id),f=(v,n=3)=>(+v.toFixed(n)).toString();
function steps(box,list){box.innerHTML="";list.forEach(s=>{if(s.h){const h=document.createElement("h3");h.textContent=s.h;box.append(h)}else{const e=document.createElement("div");e.className="eq";box.append(e);katex.render(s,e,{displayMode:true,throwOnError:false})}})}

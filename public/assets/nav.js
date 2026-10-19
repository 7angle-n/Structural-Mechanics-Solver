/* One shared top menu + footer for every page. Edit the list below and all pages update. */
(function(){
 var s=document.currentScript,r=s.getAttribute("data-root")||"",p=location.pathname.replace(/\\/g,"/");
 var L=[["Solvers","solvers/index.html","solvers"],["Beams","beam/index.html","beam"],["Stress","stress/index.html","stress"],["Trusses","truss/index.html","truss"],["Columns","columns/index.html","columns"],["Axial","axial/index.html","axial"],["Torsion","torsion/index.html","torsion"],["Inertia","inertia/index.html","inertia"],["Mohr's circle","mohr/index.html","mohr"],["Capabilities","index.html#cap",""],["Methods","index.html#methods",""],["Docs","docs/index.html","docs"]];
 var n=document.createElement("nav");n.className="nav";n.setAttribute("aria-label","Main");
 n.innerHTML='<a class="brand" href="'+r+'index.html">Structural Mechanics Solver</a>'+L.map(function(x){return '<a href="'+r+x[1]+'"'+(x[2]&&p.indexOf("/"+x[2]+"/")>-1?' aria-current="page"':"")+'>'+x[0]+'</a>'}).join("");
 s.parentNode.insertBefore(n,s);
 var m=document.createElement("meta");m.name="author";m.content="Nafis Al Zahid";document.head.appendChild(m);
 document.addEventListener("DOMContentLoaded",function(){if(!document.querySelector("footer")){var f=document.createElement("footer");f.className="sitefoot";f.textContent="Structural Mechanics Solver · created by Nafis Al Zahid";document.body.appendChild(f)}});
})();

const gamma = 1.4;
function thetaBetaEq(beta, M1, thetaRad){const s=Math.sin(beta);const c=Math.cos(beta);
const num=M1*M1*s*s-1;const den=M1*M1*(gamma+Math.cos(2*beta))+2;
return 2*(1/Math.tan(beta))*(num/den)-Math.tan(thetaRad);}
function solveBetaWeak(M1,thetaDeg){const t=thetaDeg*Math.PI/180;let lo=t+1e-6,hi=Math.PI/2-1e-6;
let flo=thetaBetaEq(lo,M1,t),fhi=thetaBetaEq(hi,M1,t);if(flo*fhi>0)return null;
for(let i=0;i<100;i++){let mid=0.5*(lo+hi);let fmid=thetaBetaEq(mid,M1,t);if(flo*fmid<=0){hi=mid;fhi=fmid;}
else{lo=mid;flo=fmid;}}return 0.5*(lo+hi);}
function obliqueShockRelations(M1,beta,theta){const Mn1=M1*Math.sin(beta);
const p2p1=1+2*gamma/(gamma+1)*(Mn1*Mn1-1);
const rho2rho1=(gamma+1)*Mn1*Mn1/((gamma-1)*Mn1*Mn1+2);
const T2T1=p2p1/rho2rho1;const Mn2=Math.sqrt((1+0.5*(gamma-1)*Mn1*Mn1)/(gamma*Mn1*Mn1-0.5*(gamma-1)));
const M2=Mn2/Math.sin(beta-theta);return{p2p1,rho2rho1,T2T1,M2};}
function arcPath(cx,cy,r,a1,a2){const x1=cx+r*Math.cos(a1),y1=cy-r*Math.sin(a1);
const x2=cx+r*Math.cos(a2),y2=cy-r*Math.sin(a2);const largeArcFlag=Math.abs(a2-a1)>Math.PI?1:0;
const sweepFlag=a2>a1?0:1;return`M ${x1.toFixed(1)},${y1.toFixed(1)} A ${r},${r} 0 ${largeArcFlag} ${sweepFlag} ${x2.toFixed(1)},${y2.toFixed(1)}`;}
function updateDiagram(thetaDeg,betaDeg){const wedge=document.getElementById("wedgeRamp"),
shock=document.getElementById("shockLine"),tArc=document.getElementById("thetaArc"),bArc=document.getElementById("betaArc"),
tLbl=document.getElementById("thetaLabel"),bLbl=document.getElementById("betaLabel");
const x0=140,y0=160,th=thetaDeg*Math.PI/180,b=betaDeg*Math.PI/180;
const Lw=80,xW=x0+Lw*Math.cos(th),yW=y0-Lw*Math.sin(th);wedge.setAttribute("x2",xW);wedge.setAttribute("y2",yW);
const Ls=130,xS=x0+Ls*Math.cos(b),yS=y0-Ls*Math.sin(b);shock.setAttribute("x2",xS);shock.setAttribute("y2",yS);
const rT=28;tArc.setAttribute("d",arcPath(x0,y0,rT,0,th));tLbl.setAttribute("x",(x0+(rT+12)*Math.cos(th/2)).toFixed(1));
tLbl.setAttribute("y",(y0-(rT+12)*Math.sin(th/2)).toFixed(1));tLbl.textContent=`θ=${thetaDeg.toFixed(1)}°`;
const rB=44;bArc.setAttribute("d",arcPath(x0,y0,rB,0,b));bLbl.setAttribute("x",(x0+(rB+12)*Math.cos(b/2)).toFixed(1));
bLbl.setAttribute("y",(y0-(rB+12)*Math.sin(b/2)).toFixed(1));bLbl.textContent=`β=${betaDeg.toFixed(1)}°`;}
function calculate(){const M1=parseFloat(document.getElementById("M1").value);
const thetaDeg=parseFloat(document.getElementById("theta").value);
const warn=document.getElementById("warning");if(isNaN(M1)||isNaN(thetaDeg)||M1<=1){
warn.textContent="M₁ must be >1 for attached shock.";updateDiagram(thetaDeg||0,thetaDeg||0);return;}
const beta=solveBetaWeak(M1,thetaDeg);if(beta===null){warn.textContent="No attached weak oblique shock for this case.";
updateDiagram(thetaDeg,thetaDeg);return;}warn.textContent="";
const theta=thetaDeg*Math.PI/180,betaDeg=beta*180/Math.PI;const{p2p1,rho2rho1,T2T1,M2}=obliqueShockRelations(M1,beta,theta);
document.getElementById("beta").textContent=betaDeg.toFixed(3);
document.getElementById("M2").textContent=M2.toFixed(3);
document.getElementById("p2p1").textContent=p2p1.toFixed(3);
document.getElementById("T2T1").textContent=T2T1.toFixed(3);
document.getElementById("rho2rho1").textContent=rho2rho1.toFixed(3);
updateDiagram(thetaDeg,betaDeg);}
document.getElementById("computeBtn").addEventListener("click",calculate);
calculate();
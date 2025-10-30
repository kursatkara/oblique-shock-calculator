// ----- Constants -----
const gamma = 1.4;

// θ–β–M equation for attached oblique shock
function thetaBetaEq(beta, M1, thetaRad) {
  const sinB = Math.sin(beta);
  const cosB = Math.cos(beta);
  const num = M1 ** 2 * sinB ** 2 - 1;
  const den = M1 ** 2 * (gamma + Math.cos(2 * beta)) + 2;
  const rhs = 2 * (1 / Math.tan(beta)) * (num / den);
  return rhs - Math.tan(thetaRad);
}

// Solve for weak-shock β using bisection
function solveBetaWeak(M1, thetaDeg) {
  const thetaRad = (thetaDeg * Math.PI) / 180;

  // lower and upper search limits for β
  let betaLo = thetaRad + 0.01;       // small offset above θ
  let betaHi = Math.PI / 2 - 1e-6;

  // if fLo and fHi same sign, gradually move betaLo up
  let fLo = thetaBetaEq(betaLo, M1, thetaRad);
  let fHi = thetaBetaEq(betaHi, M1, thetaRad);
  while (fLo * fHi > 0 && betaLo < betaHi - 0.001) {
    betaLo += 0.01;
    fLo = thetaBetaEq(betaLo, M1, thetaRad);
  }

  if (fLo * fHi > 0) return null; // still no sign change → detached

  // bisection
  for (let i = 0; i < 200; i++) {
    const betaMid = 0.5 * (betaLo + betaHi);
    const fMid = thetaBetaEq(betaMid, M1, thetaRad);
    if (fLo * fMid <= 0) {
      betaHi = betaMid;
      fHi = fMid;
    } else {
      betaLo = betaMid;
      fLo = fMid;
    }
  }
  return 0.5 * (betaLo + betaHi);
}

// Compute post-shock relations
function obliqueShockRelations(M1, beta, theta) {
  const Mn1 = M1 * Math.sin(beta);
  const p2p1 = 1 + (2 * gamma / (gamma + 1)) * (Mn1 ** 2 - 1);
  const rho2rho1 = ((gamma + 1) * Mn1 ** 2) / ((gamma - 1) * Mn1 ** 2 + 2);
  const T2T1 = p2p1 / rho2rho1;
  const Mn2 = Math.sqrt(
    (1 + 0.5 * (gamma - 1) * Mn1 ** 2) /
    (gamma * Mn1 ** 2 - 0.5 * (gamma - 1))
  );
  const M2 = Mn2 / Math.sin(beta - theta);
  return { p2p1, rho2rho1, T2T1, M2 };
}

// Draw SVG arc between angles
function arcPath(cx, cy, r, a1, a2) {
  const x1 = cx + r * Math.cos(a1);
  const y1 = cy - r * Math.sin(a1);
  const x2 = cx + r * Math.cos(a2);
  const y2 = cy - r * Math.sin(a2);
  const largeArcFlag = Math.abs(a2 - a1) > Math.PI ? 1 : 0;
  const sweepFlag = a2 > a1 ? 0 : 1;
  return `M ${x1.toFixed(1)},${y1.toFixed(1)} A ${r},${r} 0 ${largeArcFlag} ${sweepFlag} ${x2.toFixed(1)},${y2.toFixed(1)}`;
}

// Update wedge/shock diagram
function updateDiagram(thetaDeg, betaDeg) {
  const wedgeRamp = document.getElementById("wedgeRamp");
  const shockLine = document.getElementById("shockLine");
  const thetaArc = document.getElementById("thetaArc");
  const betaArc = document.getElementById("betaArc");
  const thetaLabel = document.getElementById("thetaLabel");
  const betaLabel = document.getElementById("betaLabel");

  const x0 = 140, y0 = 160;
  const thetaRad = (thetaDeg * Math.PI) / 180;
  const betaRad = (betaDeg * Math.PI) / 180;

  // wedge line
  const Lw = 80;
  const xW = x0 + Lw * Math.cos(thetaRad);
  const yW = y0 - Lw * Math.sin(thetaRad);
  wedgeRamp.setAttribute("x2", xW);
  wedgeRamp.setAttribute("y2", yW);

  // shock line
  const Ls = 130;
  const xS = x0 + Ls * Math.cos(betaRad);
  const yS = y0 - Ls * Math.sin(betaRad);
  shockLine.setAttribute("x2", xS);
  shockLine.setAttribute("y2", yS);

  // θ arc and label
  const rT = 28;
  thetaArc.setAttribute("d", arcPath(x0, y0, rT, 0, thetaRad));
  thetaLabel.setAttribute("x", x0 + (rT + 12) * Math.cos(thetaRad / 2));
  thetaLabel.setAttribute("y", y0 - (rT + 12) * Math.sin(thetaRad / 2));
  thetaLabel.textContent = `θ=${thetaDeg.toFixed(1)}°`;

  // β arc and label
  const rB = 44;
  betaArc.setAttribute("d", arcPath(x0, y0, rB, 0, betaRad));
  betaLabel.setAttribute("x", x0 + (rB + 12) * Math.cos(betaRad / 2));
  betaLabel.setAttribute("y", y0 - (rB + 12) * Math.sin(betaRad / 2));
  betaLabel.textContent = `β=${betaDeg.toFixed(1)}°`;
}

// Compute results and update UI
function calculate() {
  const M1 = parseFloat(document.getElementById("M1").value);
  const thetaDeg = parseFloat(document.getElementById("theta").value);
  const warning = document.getElementById("warning");

  if (isNaN(M1) || isNaN(thetaDeg) || M1 <= 1) {
    warning.textContent = "M₁ must be > 1 for an attached oblique shock.";
    return;
  }

  const beta = solveBetaWeak(M1, thetaDeg);
  if (beta === null) {
    warning.textContent = "No attached weak oblique shock for this (M₁, θ).";
    updateDiagram(thetaDeg, thetaDeg);
    return;
  }
  warning.textContent = "";

  const theta = (thetaDeg * Math.PI) / 180;
  const betaDeg = (beta * 180) / Math.PI;

  const { p2p1, rho2rho1, T2T1, M2 } = obliqueShockRelations(M1, beta, theta);

  document.getElementById("beta").textContent = betaDeg.toFixed(6);
  document.getElementById("M2").textContent = M2.toFixed(6);
  document.getElementById("p2p1").textContent = p2p1.toFixed(6);
  document.getElementById("T2T1").textContent = T2T1.toFixed(6);
  document.getElementById("rho2rho1").textContent = rho2rho1.toFixed(6);

  updateDiagram(thetaDeg, betaDeg);
}

document.getElementById("computeBtn").addEventListener("click", calculate);
calculate();


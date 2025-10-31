// ------------------------------------------------------------
// Oblique Shock Calculator
// Analytical theta–beta–Mach solution (Rudd & Lewis, 1998)
// with corrected wedge–shock geometry visualization
// ------------------------------------------------------------

const gamma = 1.4;

// Analytical θ–β–M relation (Rudd & Lewis 1998)
function betaAnalytical(M1, thetaDeg, gamma = 1.4, n = 0) {
  const theta = (thetaDeg * Math.PI) / 180;
  const mu = Math.asin(1 / M1);
  const c = Math.tan(mu) ** 2;
  const a = ((gamma - 1) / 2 + ((gamma + 1) / 2) * c) * Math.tan(theta);
  const b = ((gamma + 1) / 2 + ((gamma + 3) / 2) * c) * Math.tan(theta);
  const numerator = 4 * (1 - 3 * a * b) ** 3;
  const denominator = (27 * a * a * c + 9 * a * b - 2) ** 2;
  const d = Math.sqrt(numerator / denominator - 1);

  const term1 = (b + 9 * a * c) / (2 * (1 - 3 * a * b));
  const term2 =
    (d * (27 * a * a * c + 9 * a * b - 2)) /
    (6 * a * (1 - 3 * a * b));

  const angle = n * Math.PI / 3 + (1 / 3) * Math.atan(1 / d);
  const betaRad = Math.atan(term1 - term2 * Math.tan(angle));

  return betaRad; // radians
}

// Compute post-shock flow properties
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

// Utility: create SVG arc for labeling θ and β
function arcPath(cx, cy, r, a1, a2) {
  const x1 = cx + r * Math.cos(a1);
  const y1 = cy - r * Math.sin(a1);
  const x2 = cx + r * Math.cos(a2);
  const y2 = cy - r * Math.sin(a2);
  const largeArcFlag = Math.abs(a2 - a1) > Math.PI ? 1 : 0;
  const sweepFlag = a2 > a1 ? 0 : 1;
  return `M ${x1.toFixed(1)},${y1.toFixed(1)} A ${r},${r} 0 ${largeArcFlag} ${sweepFlag} ${x2.toFixed(1)},${y2.toFixed(1)}`;
}

// Corrected wedge–shock geometry
function updateDiagram(thetaDeg, betaDeg) {
  const wedgeBase = document.getElementById("wedgeBase");
  const wedgeRamp = document.getElementById("wedgeRamp");
  const shockLine = document.getElementById("shockLine");
  const thetaArc = document.getElementById("thetaArc");
  const betaArc = document.getElementById("betaArc");
  const thetaLabel = document.getElementById("thetaLabel");
  const betaLabel = document.getElementById("betaLabel");

  const x0 = 280;
  const y0 = 320;

  const thetaRad = (thetaDeg * Math.PI) / 180;
  const betaRad = (betaDeg * Math.PI) / 180;

  const Lbase = 200;
  const Lwedge = 120;
  const Lshock = 200;

  // Wedge baseline (horizontal)
  const xBase2 = x0 + Lbase;
  const yBase2 = y0;
  wedgeBase.setAttribute("x1", x0);
  wedgeBase.setAttribute("y1", y0);
  wedgeBase.setAttribute("x2", xBase2);
  wedgeBase.setAttribute("y2", yBase2);

  // Wedge ramp
  const xRamp2 = x0 + Lwedge * Math.cos(thetaRad);
  const yRamp2 = y0 - Lwedge * Math.sin(thetaRad);
  wedgeRamp.setAttribute("x1", x0);
  wedgeRamp.setAttribute("y1", y0);
  wedgeRamp.setAttribute("x2", xRamp2);
  wedgeRamp.setAttribute("y2", yRamp2);

  // Shock line
  const xShock2 = x0 + Lshock * Math.cos(betaRad);
  const yShock2 = y0 - Lshock * Math.sin(betaRad);
  shockLine.setAttribute("x1", x0);
  shockLine.setAttribute("y1", y0);
  shockLine.setAttribute("x2", xShock2);
  shockLine.setAttribute("y2", yShock2);

  // θ arc and label
  const rTheta = 40;
  thetaArc.setAttribute("d", arcPath(x0, y0, rTheta, 0, thetaRad));
  const thetaMid = thetaRad / 2;
  thetaLabel.setAttribute("x", x0 + (rTheta + 20) * Math.cos(thetaMid));
  thetaLabel.setAttribute("y", y0 - (rTheta + 20) * Math.sin(thetaMid));
  thetaLabel.textContent = `θ=${thetaDeg.toFixed(1)}°`;

  // β arc and label
  const rBeta = 60;
  betaArc.setAttribute("d", arcPath(x0, y0, rBeta, 0, betaRad));
  const betaMid = betaRad / 2;
  betaLabel.setAttribute("x", x0 + (rBeta + 20) * Math.cos(betaMid));
  betaLabel.setAttribute("y", y0 - (rBeta + 20) * Math.sin(betaMid));
  betaLabel.textContent = `β=${betaDeg.toFixed(1)}°`;
}

// Main computation
function calculate() {
  const M1 = parseFloat(document.getElementById("M1").value);
  const thetaDeg = parseFloat(document.getElementById("theta").value);
  const warning = document.getElementById("warning");

  if (isNaN(M1) || isNaN(thetaDeg) || M1 <= 1) {
    warning.textContent = "M₁ must be > 1 for an attached oblique shock.";
    return;
  }

  const betaRad = betaAnalytical(M1, thetaDeg, gamma, 0);
  const betaDeg = betaRad * 180 / Math.PI;

  if (isNaN(betaRad) || betaDeg <= thetaDeg || betaDeg >= 90) {
    warning.textContent = "No attached weak oblique shock for this (M₁, θ).";
    updateDiagram(thetaDeg, thetaDeg);
    return;
  }

  warning.textContent = "";

  const thetaRad = thetaDeg * Math.PI / 180;
  const { p2p1, rho2rho1, T2T1, M2 } = obliqueShockRelations(M1, betaRad, thetaRad);

  // Output with precision
  document.getElementById("beta").textContent = betaDeg.toFixed(2);
  document.getElementById("M2").textContent = M2.toFixed(4);
  document.getElementById("p2p1").textContent = p2p1.toFixed(4);
  document.getElementById("T2T1").textContent = T2T1.toFixed(4);
  document.getElementById("rho2rho1").textContent = rho2rho1.toFixed(4);

  updateDiagram(thetaDeg, betaDeg);
}

document.getElementById("computeBtn").addEventListener("click", calculate);
calculate();


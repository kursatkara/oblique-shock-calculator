const gamma = 1.4;

// θ–β–M relation (Rudd & Lewis)
function betaAnalytical(M1, thetaDeg, gamma = 1.4, n = 0) {
  const theta = (thetaDeg * Math.PI) / 180;
  const mu = Math.asin(1 / M1);
  const c = Math.tan(mu) ** 2;
  const a = ((gamma - 1) / 2 + ((gamma + 1) / 2) * c) * Math.tan(theta);
  const b = ((gamma + 1) / 2 + ((gamma + 3) / 2) * c) * Math.tan(theta);
  const num = 4 * (1 - 3 * a * b) ** 3;
  const den = (27 * a * a * c + 9 * a * b - 2) ** 2;
  const d = Math.sqrt(num / den - 1);
  const term1 = (b + 9 * a * c) / (2 * (1 - 3 * a * b));
  const term2 = (d * (27 * a * a * c + 9 * a * b - 2)) / (6 * a * (1 - 3 * a * b));
  const angle = (n * Math.PI) / 3 + (1 / 3) * Math.atan(1 / d);
  return Math.atan(term1 - term2 * Math.tan(angle));
}

function obliqueShock(M1, beta, theta) {
  const Mn1 = M1 * Math.sin(beta);
  const p2p1 = 1 + (2 * gamma / (gamma + 1)) * (Mn1 ** 2 - 1);
  const rho2rho1 = ((gamma + 1) * Mn1 ** 2) / ((gamma - 1) * Mn1 ** 2 + 2);
  const T2T1 = p2p1 / rho2rho1;
  const Mn2 = Math.sqrt((1 + 0.5 * (gamma - 1) * Mn1 ** 2) / (gamma * Mn1 ** 2 - 0.5 * (gamma - 1)));
  const M2 = Mn2 / Math.sin(beta - theta);
  const p02p01 = p2p1 * ((1 + 0.5 * (gamma - 1) * Mn2 ** 2) / (1 + 0.5 * (gamma - 1) * Mn1 ** 2)) ** (gamma / (gamma - 1));
  return { Mn1, Mn2, M2, p2p1, rho2rho1, T2T1, p02p01 };
}

function arcPath(cx, cy, r, a1, a2) {
  const x1 = cx + r * Math.cos(a1);
  const y1 = cy - r * Math.sin(a1);
  const x2 = cx + r * Math.cos(a2);
  const y2 = cy - r * Math.sin(a2);
  return `M ${x1},${y1} A ${r},${r} 0 0 0 ${x2},${y2}`;
}

function updateDiagram(thetaDeg, betaDeg) {
  const wedgeBase = document.getElementById("wedgeBase");
  const wedgeRamp = document.getElementById("wedgeRamp");
  const shockLine = document.getElementById("shockLine");
  const thetaArc = document.getElementById("thetaArc");
  const betaArc = document.getElementById("betaArc");
  const thetaLabel = document.getElementById("thetaLabel");
  const betaLabel = document.getElementById("betaLabel");

  const x0 = 280, y0 = 320;
  const thetaRad = (thetaDeg * Math.PI) / 180;
  const betaRad = (betaDeg * Math.PI) / 180;

  const Lbase = 220, Lwedge = 150, Lshock = 220;

  wedgeBase.setAttribute("x1", x0);
  wedgeBase.setAttribute("y1", y0);
  wedgeBase.setAttribute("x2", x0 + Lbase);
  wedgeBase.setAttribute("y2", y0);

  wedgeRamp.setAttribute("x1", x0);
  wedgeRamp.setAttribute("y1", y0);
  wedgeRamp.setAttribute("x2", x0 + Lwedge * Math.cos(thetaRad));
  wedgeRamp.setAttribute("y2", y0 - Lwedge * Math.sin(thetaRad));

  shockLine.setAttribute("x1", x0);
  shockLine.setAttribute("y1", y0);
  shockLine.setAttribute("x2", x0 + Lshock * Math.cos(betaRad));
  shockLine.setAttribute("y2", y0 - Lshock * Math.sin(betaRad));

  thetaArc.setAttribute("d", arcPath(x0, y0, 40, 0, thetaRad));
  betaArc.setAttribute("d", arcPath(x0, y0, 60, 0, betaRad));

  thetaLabel.setAttribute("x", x0 + 60 * Math.cos(thetaRad / 2));
  thetaLabel.setAttribute("y", y0 - 60 * Math.sin(thetaRad / 2));
  thetaLabel.textContent = `θ=${thetaDeg.toFixed(1)}°`;

  betaLabel.setAttribute("x", x0 + 85 * Math.cos(betaRad / 2));
  betaLabel.setAttribute("y", y0 - 85 * Math.sin(betaRad / 2));
  betaLabel.textContent = `β=${betaDeg.toFixed(1)}°`;
}

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
  const { Mn1, Mn2, M2, p2p1, rho2rho1, T2T1, p02p01 } = obliqueShock(M1, betaRad, thetaRad);

  document.getElementById("thetaOut").textContent = thetaDeg.toFixed(1);
  document.getElementById("beta").textContent = betaDeg.toFixed(2);
  document.getElementById("M1n").textContent = Mn1.toFixed(4);
  document.getElementById("M2n").textContent = Mn2.toFixed(4);
  document.getElementById("M2").textContent = M2.toFixed(4);
  document.getElementById("p2p1").textContent = p2p1.toFixed(4);
  document.getElementById("rho2rho1").textContent = rho2rho1.toFixed(4);
  document.getElementById("T2T1").textContent = T2T1.toFixed(4);
  document.getElementById("p02p01").textContent = p02p01.toFixed(4);

  updateDiagram(thetaDeg, betaDeg);
}

document.getElementById("computeBtn").addEventListener("click", calculate);
calculate();

// ------------------------------------------------------------
// Oblique Shock Calculator
// Uses analytical θ–β–M solution (Rudd & Lewis, AIAA J. Aircraft, 1998)
// Computes full set of shock properties for γ=1.4
// ------------------------------------------------------------

const gamma = 1.4;

// 1. Analytical θ–β–M relation
// n = 0 → weak shock, n = 1 → strong shock
function betaAnalytical(M1, thetaDeg, gamma = 1.4, n = 0) {
  const theta = (thetaDeg * Math.PI) / 180;
  const mu = Math.asin(1 / M1); // Mach angle
  const c = Math.tan(mu) ** 2;

  const a =
    ((gamma - 1) / 2 + ((gamma + 1) / 2) * c) * Math.tan(theta);
  const b =
    ((gamma + 1) / 2 + ((gamma + 3) / 2) * c) * Math.tan(theta);

  const numerator = 4 * (1 - 3 * a * b) ** 3;
  const denominator =
    (27 * a * a * c + 9 * a * b - 2) ** 2;
  const d = Math.sqrt(numerator / denominator - 1);

  const term1 =
    (b + 9 * a * c) / (2 * (1 - 3 * a * b));
  const term2 =
    (d * (27 * a * a * c + 9 * a * b - 2)) /
    (6 * a * (1 - 3 * a * b));

  const angle =
    (n * Math.PI) / 3 + (1 / 3) * Math.atan(1 / d);

  const betaRad = Math.atan(
    term1 - term2 * Math.tan(angle)
  );

  return betaRad; // radians
}

// 2. Normal-shock property relations (component normal to the shock)
function obliqueShockRelations(M1, beta, theta) {
  // Upstream Mach normal to the shock
  const Mn1 = M1 * Math.sin(beta);

  // Static pressure ratio across the shock
  const p2p1 =
    1 +
    (2 * gamma / (gamma + 1)) * (Mn1 ** 2 - 1);

  // Density ratio
  const rho2rho1 =
    ((gamma + 1) * Mn1 ** 2) /
    ((gamma - 1) * Mn1 ** 2 + 2);

  // Temperature ratio
  const T2T1 = p2p1 / rho2rho1;

  // Downstream normal Mach
  const Mn2 = Math.sqrt(
    (1 +
      0.5 * (gamma - 1) * Mn1 ** 2) /
      (gamma * Mn1 ** 2 -
        0.5 * (gamma - 1))
  );

  // Downstream Mach in the flow direction after turning by θ
  const M2 = Mn2 / Math.sin(beta - theta);

  // Stagnation pressure ratio across shock (p02/p01)
  // p0 = p * (1 + (γ-1)/2 M^2)^{γ/(γ-1)}
  // For a normal shock, ratio is:
  // p02/p01 = [p2/p1] * [(1 + (γ-1)/2 * Mn2^2) / (1 + (γ-1)/2 * Mn1^2)]^{γ/(γ-1)}
  const termUp =
    1 + ((gamma - 1) / 2) * Mn1 ** 2;
  const termDown =
    1 + ((gamma - 1) / 2) * Mn2 ** 2;
  const p02p01 =
    p2p1 *
    (termDown / termUp) ** (gamma / (gamma - 1));

  return {
    Mn1,
    Mn2,
    M2,
    p2p1,
    rho2rho1,
    T2T1,
    p02p01
  };
}

// 3. SVG arc path helper for θ and β
function arcPath(cx, cy, r, a1, a2) {
  const x1 = cx + r * Math.cos(a1);
  const y1 = cy - r * Math.sin(a1);
  const x2 = cx + r * Math.cos(a2);
  const y2 = cy - r * Math.sin(a2);
  const largeArcFlag = Math.abs(a2 - a1) > Math.PI ? 1 : 0;
  const sweepFlag = a2 > a1 ? 0 : 1;
  return `M ${x1.toFixed(1)},${y1.toFixed(
    1
  )} A ${r},${r} 0 ${largeArcFlag} ${sweepFlag} ${x2.toFixed(
    1
  )},${y2.toFixed(1)}`;
}

// 4. Draw wedge, shock, and labels
function updateDiagram(thetaDeg, betaDeg) {
  const wedgeBase = document.getElementById("wedgeBase");
  const wedgeRamp = document.getElementById("wedgeRamp");
  const shockLine = document.getElementById("shockLine");
  const thetaArc = document.getElementById("thetaArc");
  const betaArc = document.getElementById("betaArc");
  const thetaLabel = document.getElementById("thetaLabel");
  const betaLabel = document.getElementById("betaLabel");

  // Anchor (wedge tip / shock origin)
  const x0 = 280;
  const y0 = 320;

  const thetaRad = (thetaDeg * Math.PI) / 180;
  const betaRad = (betaDeg * Math.PI) / 180;

  // Drawing lengths
  const Lbase = 220;
  const Lwedge = 150;
  const Lshock = 220;

  // Baseline along freestream / lower wedge surface
  const xBase2 = x0 + Lbase;
  const yBase2 = y0;
  wedgeBase.setAttribute("x1", x0);
  wedgeBase.setAttribute("y1", y0);
  wedgeBase.setAttribute("x2", xBase2);
  wedgeBase.setAttribute("y2", yBase2);

  // Upper wedge face deflected by θ
  const xRamp2 = x0 + Lwedge * Math.cos(thetaRad);
  const yRamp2 = y0 - Lwedge * Math.sin(thetaRad);
  wedgeRamp.setAttribute("x1", x0);
  wedgeRamp.setAttribute("y1", y0);
  wedgeRamp.setAttribute("x2", xRamp2);
  wedgeRamp.setAttribute("y2", yRamp2);

  // Shock at β
  const xShock2 = x0 + Lshock * Math.cos(betaRad);
  const yShock2 = y0 - Lshock * Math.sin(betaRad);
  shockLine.setAttribute("x1", x0);
  shockLine.setAttribute("y1", y0);
  shockLine.setAttribute("x2", xShock2);
  shockLine.setAttribute("y2", yShock2);

  // θ arc (between baseline and wedge face)
  const rTheta = 40;
  thetaArc.setAttribute(
    "d",
    arcPath(x0, y0, rTheta, 0, thetaRad)
  );

  // β arc (between baseline and shock)
  const rBeta = 60;
  betaArc.setAttribute(
    "d",
    arcPath(x0, y0, rBeta, 0, betaRad)
  );

  // θ label
  const thetaMid = thetaRad / 2;
  const thetaLabelX =
    x0 + (rTheta + 22) * Math.cos(thetaMid);
  const thetaLabelY =
    y0 - (rTheta + 22) * Math.sin(thetaMid);
  thetaLabel.setAttribute("x", thetaLabelX.toFixed(1));
  thetaLabel.setAttribute("y", thetaLabelY.toFixed(1));
  thetaLabel.setAttribute("text-anchor", "middle");
  thetaLabel.setAttribute(
    "dominant-baseline",
    "middle"
  );
  thetaLabel.textContent =
    "θ=" + thetaDeg.toFixed(1) + "°";

  // β label
  const betaMid = betaRad / 2;
  const betaLabelX =
    x0 + (rBeta + 22) * Math.cos(betaMid);
  const betaLabelY =
    y0 - (rBeta + 22) * Math.sin(betaMid);
  betaLabel.setAttribute("x", betaLabelX.toFixed(1));
  betaLabel.setAttribute("y", betaLabelY.toFixed(1));
  betaLabel.setAttribute("text-anchor", "middle");
  betaLabel.setAttribute(
    "dominant-baseline",
    "middle"
  );
  betaLabel.textContent =
    "β=" + betaDeg.toFixed(1) + "°";
}

// 5. Main calculation and DOM update
function calculate() {
  const M1 = parseFloat(document.getElementById("M1").value);
  const thetaDeg = parseFloat(document.getElementById("theta").value);
  const warning = document.getElementById("warning");

  // sanity check
  if (isNaN(M1) || isNaN(thetaDeg) || M1 <= 1) {
    warning.textContent =
      "M₁ must be > 1 for an attached oblique shock.";
    return;
  }

  // Weak-shock branch
  const betaRad = betaAnalytical(M1, thetaDeg, gamma, 0);
  const betaDeg = (betaRad * 180) / Math.PI;

  // Physical validity check
  if (
    isNaN(betaRad) ||
    betaDeg <= thetaDeg ||
    betaDeg >= 90
  ) {
    warning.textContent =
      "No attached weak oblique shock for this (M₁, θ).";
    updateDiagram(thetaDeg, thetaDeg); // draw something consistent
    // Clear outputs
    document.getElementById("thetaOut").textContent = "–";
    document.getElementById("beta").textContent = "–";
    document.getElementById("M1n").textContent = "–";
    document.getElementById("M2n").textContent = "–";
    document.getElementById("M2").textContent = "–";
    document.getElementById("p2p1").textContent = "–";
    document.getElementById("rho2rho1").textContent = "–";
    document.getElementById("T2T1").textContent = "–";
    document.getElementById("p02p01").textContent = "–";
    return;
  }

  warning.textContent = "";

  const thetaRad = (thetaDeg * Math.PI) / 180;
  const {
    Mn1,
    Mn2,
    M2,
    p2p1,
    rho2rho1,
    T2T1,
    p02p01
  } = obliqueShockRelations(M1, betaRad, thetaRad);

  // Update numeric results
  document.getElementById("thetaOut").textContent =
    thetaDeg.toFixed(1);
  document.getElementById("beta").textContent =
    betaDeg.toFixed(2);
  document.getElementById("M1n").textContent =
    Mn1.toFixed(4);
  document.getElementById("M2n").textContent =
    Mn2.toFixed(4);
  document.getElementById("M2").textContent =
    M2.toFixed(4);
  document.getElementById("p2p1").textContent =
    p2p1.toFixed(4);
  document.getElementById("rho2rho1").textContent =
    rho2rho1.toFixed(4);
  document.getElementById("T2T1").textContent =
    T2T1.toFixed(4);
  document.getElementById("p02p01").textContent =
    p02p01.toFixed(4);

  // Update diagram
  updateDiagram(thetaDeg, betaDeg);
}

// Hook up button + initial run
document
  .getElementById("computeBtn")
  .addEventListener("click", calculate);

calculate();

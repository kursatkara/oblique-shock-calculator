// Constants
const gamma = 1.4;

// Compute oblique shock relations
function computeObliqueShock(M1, thetaDeg) {
  const theta = thetaDeg * Math.PI / 180;
  
  // Iteratively solve for beta using Newton-Raphson
  let beta = theta + Math.PI / 10; // initial guess
  for (let i = 0; i < 100; i++) {
    const f = 2 * Math.tan(beta - theta) *
      ((M1 * M1 * Math.sin(beta) ** 2 - 1) /
      (M1 * M1 * (gamma + Math.cos(2 * beta)) + 2)) - Math.tan(theta);
    const df = (2 * (M1 * M1 * Math.sin(beta) * Math.cos(beta) * (gamma + Math.cos(2 * beta)) -
      (M1 * M1 * Math.sin(beta) ** 2 - 1) * (2 * M1 * M1 * Math.sin(beta) * Math.cos(beta))) /
      ((M1 * M1 * (gamma + Math.cos(2 * beta)) + 2) ** 2));
    beta -= f / df;
  }

  const Mn1 = M1 * Math.sin(beta);
  const Mn2 = Math.sqrt((1 + (gamma - 1) / 2 * Mn1 ** 2) / (gamma * Mn1 ** 2 - (gamma - 1) / 2));
  const M2 = Mn2 / Math.sin(beta - theta);
  const p2p1 = 1 + 2 * gamma / (gamma + 1) * (Mn1 ** 2 - 1);
  const T2T1 = p2p1 * ((gamma + 1) * Mn1 ** 2) / (2 + (gamma - 1) * Mn1 ** 2);
  const rho2rho1 = p2p1 / T2T1;

  return { beta: beta * 180 / Math.PI, M2, p2p1, T2T1, rho2rho1 };
}

function calculate() {
  const M1 = parseFloat(document.getElementById("M1").value);
  const theta = parseFloat(document.getElementById("theta").value);
  if (isNaN(M1) || isNaN(theta)) return;

  const res = computeObliqueShock(M1, theta);
  document.getElementById("beta").textContent = res.beta.toFixed(3);
  document.getElementById("M2").textContent = res.M2.toFixed(3);
  document.getElementById("p2p1").textContent = res.p2p1.toFixed(3);
  document.getElementById("T2T1").textContent = res.T2T1.toFixed(3);
  document.getElementById("rho2rho1").textContent = res.rho2rho1.toFixed(3);
}


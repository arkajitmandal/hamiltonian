'use strict';
//const math = require('mathjs');
//const numeric = require('../numeric.js');

// Class for creating a Basis Function
class BasisFunction {
  constructor(expression, up, down) {
    up = up || 1.0;
    down = down || -1.0;
    this.fx = expression;
    this.upperLimit = up;
    this.lowerLimit = down;
  }
  // Set Limit
  setLimits(limit) {
    this.upperLimit = math.max(limit);
    this.lowerLimit = math.min(limit);
  }
  // Evaluate the function at a given x
  at(x, fx = this.fx) {
    return math.eval(fx, { x: x });
  }
  // Derivative arround x
  derivative(x, analytical = false) {
    // Numerical Derivative
    if (!analytical) {
      let h = 0.0000001;
      return (this.at(x + h) - this.at(x - h)) / (2 * h);
    }
    // Analytical Derivative
    let Dx = math.eval(math.derivative(this.fx, 'x').toString(), { x: x });
    return Dx;
  }
  // Integrate the function
  // < psi | psi >
  integrate(from = this.lowerLimit, to = this.upperLimit, points = 100) {
    let dx = (to - from) / points;
    let total = 0.0;
    for (var i = 0; i < points; i++) {
      let x = from + i * dx;
      let fx2 = '(' + this.fx + ')^2';
      total += dx / 2 * (this.at(x, fx2) + this.at(x + dx, fx2));
    }
    return total;
  }
  // Normalize the function
  // also return the normalization constant
  //         ______________
  // 1/N =  V < psi | psi >
  normalize(from = this.lowerLimit, to = this.upperLimit) {
    let N = math.sqrt(this.integrate(from, to));
    this.fx = (1 / N).toString() + '* (' + this.fx + ')';
    this.upperLimit = to;
    this.lowerLimit = from;
    return N;
  }
  // Double Derivative
  // d^2/dx^2 (psi)
  doubleDerivative(x, analytical = false) {
    // Numerical 5 point Derivative
    if (!analytical) {
      let h = 0.000001;
      let D1 = this.derivative(x + 2 * h);
      let D2 = this.derivative(x + h);
      let D3 = this.derivative(x - h);
      let D4 = this.derivative(x - 2 * h);
      return (D4 - D1 + 8 * (D2 - D3)) / (12 * h);
    }
    let Dx = math.derivative(this.fx, 'x').toString();
    let D2x = math.derivative(Dx, 'x').toString();
    return math.eval(D2x, { x: x });
  }
  // Kinetic Energy
  // - < psi | d^2/dx^2 | psi >/m
  kinetic(m = 1.0, analytical = false, points = 100) {
    // Value of Kinetic Energy
    let dx = (this.upperLimit - this.lowerLimit) / points;
    let total = 0.0;
    for (var i = 0; i < points; i++) {
      let x = this.lowerLimit + i * dx;
      let D2x0 = this.doubleDerivative(x, analytical);
      let D2x1 = this.doubleDerivative(x + dx, analytical);
      let fx0 = this.at(x);
      let fx1 = this.at(x + dx);
      total += dx / 2 * (D2x0 * fx0 + D2x1 * fx1);
    }
    return -total / m;
  }
}

// Class for creating Basis Set
class Basis {
  constructor(basisSet) {
    this.set = basisSet;
    this.upperLimit = basisSet[0].upperLimit;
    this.lowerLimit = basisSet[0].lowerLimit;
    this.normalized = false;
  }
  // Normalize all functions
  normalize(lower = this.lowerLimit, upper = this.upperLimit) {
    for (var i = 0; i < this.set.length; i++) {
      this.set[i].normalize(lower, upper);
    }
    this.upperLimit = upper;
    this.lowerLimit = lower;
    this.normalized = true;
  }
  // Overlap between ith and jth psi
  // instead of pass index, psi itself
  // can be passed through func = [psi1,psi2]
  overlap(i, j, func = false, points = 100) {
    let psi1;
    let psi2;
    if (func === false) {
      psi1 = this.set[i];
      psi2 = this.set[j];
    } else {
      psi1 = func[0];
      psi2 = func[1];
    }
    let dx = (this.upperLimit - this.lowerLimit) / points;
    let total = 0.0;
    for (var k = 0; k < points; k++) {
      let x = this.lowerLimit + k * dx;
      total += dx / 2 * (psi1.at(x) * psi2.at(x) + psi1.at(x + dx) * psi2.at(x + dx));
    }
    return total;
  }
  // 2nd order derivative coupling
  //  < i | d^2/dx^2 | j >
  Dij(i, j, analytical = false, points = 100) {
    let dx = (this.upperLimit - this.lowerLimit) / points;
    let total = 0.0;
    for (var k = 0; k < points; k++) {
      let x = this.lowerLimit + k * dx;
      let D2x0 = this.set[j].doubleDerivative(x, analytical);
      let D2x1 = this.set[j].doubleDerivative(x + dx, analytical);
      let fx0 = this.set[i].at(x);
      let fx1 = this.set[i].at(x + dx);
      total += dx / 2 * (D2x0 * fx0 + D2x1 * fx1);
    }
    return total;
  }
  // Gram-Smidth Orthonormalization
  orthoNormalize() {
    this.set[0].normalize();
    let orthoBasis = [this.set[0]];
    let _;
    for (var i = 1; i < this.set.length; i++) {
      // Get ith function
      let thisFunction = '(' + this.set[i].fx + ')';
      for (var j = i - 1; j >= 0; j--) {
        let Sij = this.overlap(_, _, [this.set[i], orthoBasis[j]]);
        let Sjj = this.overlap(_, _, [orthoBasis[j], orthoBasis[j]]);
        thisFunction += '- (' + (Sij / Sjj).toString() + ') * (' + orthoBasis[j].fx + ')';
      }
      orthoBasis.push(new BasisFunction(thisFunction, this.upperLimit, this.lowerLimit));
      orthoBasis[i].normalize();
    }
    this.set = orthoBasis;
  }
}

// Evaluate Hamiltonian Matrix Elements
// V is the state indepent potential
// By default it assumes non-ortho-normalized
// basis, should use orthgonal = true
// if the basis is ortho-normal for faster evaluation
function hamiltonian(basis, m = 1, V = 0.0, orthgonal = false) {
  // Empty NxN array for total Hamiltonian
  let N = basis.set.length;
  let H = Array.from({ length: N }, () => new Array(N).fill(0));
  for (var i = 0; i < N; i++) {
    for (var j = 0; j < N; j++) {
      // eslint-disable-next-line new-cap
      H[i][j] = -basis.Dij(i, j) / m;
      if (orthgonal) {
        H[i][j] += V * basis.overlap(i, j);
      } else if (i === j) {
        H[i][i] += V;
      }
    }
  }
  return H;
}

// Evaluate eigenvalue
// and eigen vector
// Assuming they are real
function eigen(H) {
  let eig = numeric.eig(H);
  return [eig.lambda.x, eig.E.x];
}

// Generate some common Basis Function
function commonBasis(name, N, up = 1, down = -1) {
  if (name === 'legendre') {
    let fx = [];
    let bfx = Array(N).fill('');
    for (var n = 0; n < N; n++) {
      bfx[n] = '0';
      // The ith polynomial generated here
      for (var k = 0; k <= n; k++) {
        let coeff = math.combinations(n, k) * math.combinations(n + k, k);
        bfx[n] =
          ' (' + coeff.toString() + ' * ( x/2 - 1/2 )^' + k.toString() + ' ) + ' + bfx[n];
      }
      fx.push(new BasisFunction(bfx[n], up, down));
    }
    let legendre = new Basis(fx);
    return legendre;
  }
}

module.exports = {
  BasisFunction: BasisFunction,
  Basis: Basis,
  hamiltonian: hamiltonian,
  eigen: eigen,
  commonBasis: commonBasis
};

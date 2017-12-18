const assert = require('assert');
const hamiltonian = require('../index.js');
const math = require('mathjs');

describe('Basis Function', () => {
  it('Creation of a Basis Function', () => {
    let fx = new hamiltonian.BasisFunction('x + 2');
    assert(fx.at(9.0) === 11.0, 'function not created properly');
  });
  it('Derivative of a Basis Function', () => {
    let fx = new hamiltonian.BasisFunction('x^4');
    assert(fx.derivative(0.0, false) === 0.0, 'derivative not 0.0');
    assert(
      fx.derivative(0.0, false) === fx.derivative(0.0, true),
      'analytical and numerical should match at this point'
    );
    assert(
      Math.abs(fx.derivative(5.0, false) - fx.derivative(5.0, true)) < 0.00001,
      ' two derivatives should be close'
    );
  });
  it('Integration of a Basis Function', () => {
    let fx;
    fx = new hamiltonian.BasisFunction('x^2');
    let intfx = math.eval('x^5/5', { x: 2.0 }) - math.eval('x^5/5', { x: -2.0 });
    assert(
      Math.abs(fx.integrate(-2.0, 2.0) / intfx) > 0.9999,
      'Integration is far from analytical value'
    );
    assert(
      Math.abs(fx.integrate(-2.0, 2.0) / intfx) < 1.0001,
      'Integration is far from analytical value'
    );
  });
  it('Normalization of a Basis Function', () => {
    let fx;
    fx = new hamiltonian.BasisFunction('x^2 + x^7');
    fx.normalize(-5.0, 5.0);
    assert(Math.abs(fx.integrate(-5.0, 5.0) - 1.0) < 0.00001, 'Function not normalized');
    fx = new hamiltonian.BasisFunction('sin(pi * x/2.0)');
    let N = fx.normalize(-1.0, 1.0);
    assert(Math.abs(N - 1.0) < 0.00001, 'Particle in  a box Function not normalized');
  });
  it('Double Derivative of a Basis Function', () => {
    let fx;
    fx = new hamiltonian.BasisFunction('x^2');
    let D2x = fx.doubleDerivative(0.0);
    assert(Math.abs(D2x - 2.0) < 0.001, 'Double Derivative should be zero');
    fx = new hamiltonian.BasisFunction('sin(x)');
    D2x = fx.doubleDerivative(2.0);
    assert(
      Math.abs(D2x - fx.doubleDerivative(2.0, true)) < 0.01,
      'analytical and numerical are far'
    );
  });
  it('Kinetic Energy operator', () => {
    let fx;
    fx = new hamiltonian.BasisFunction('sin(x)');
    fx.normalize();
    let K = fx.kinetic();
    assert(Math.abs(K - 1.0) < 0.001, 'Kinetic energy should be 1');
    K = fx.kinetic(undefined, true);
    assert(Math.abs(K - 1.0) < 0.001, 'Kinetic energy should be 1');
  });
  it('Set Limit', () => {
    let fx;
    fx = new hamiltonian.BasisFunction('sin(x)');
    fx.setLimits([-2.0, -5.0]);
    assert(fx.lowerLimit === -5.0, 'wrong lower limits set');
    assert(fx.upperLimit === -2.0, 'wrong upper limits set');
  });
});
describe('Basis Set', () => {
  it('Creation of a Basis Set', () => {
    let fx1 = new hamiltonian.BasisFunction('x + 2');
    let fx2 = new hamiltonian.BasisFunction('sin(x) + 2');
    let basis = new hamiltonian.Basis([fx1, fx2]);
    assert(basis.set[0].fx === fx1.fx, 'Basis not created properly');
    assert(basis.set[1].fx === fx2.fx, 'Basis not created properly');
  });
  it('Overlap between functions', () => {
    let fx1 = new hamiltonian.BasisFunction('x + 2');
    let fx2 = new hamiltonian.BasisFunction('sin(x) + 2');
    let basis = new hamiltonian.Basis([fx1, fx2]);
    assert(basis.overlap(0, 1) === basis.overlap(1, 0), 'Overlap must be symmetric');
    assert(basis.overlap(0, 1) !== 0.0, 'Overlap should not be 0');
    fx1 = new hamiltonian.BasisFunction('sin(x)');
    fx2 = new hamiltonian.BasisFunction('cos(x)');
    basis = new hamiltonian.Basis([fx1, fx2]);
    assert(basis.overlap(0, 1) <= 0.000001, 'Overlap should be 0');
  });
});

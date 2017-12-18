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
      Math.abs(fx.integrate(-2.0, 2.0) / intfx) > 0.999,
      'Integration is far from analytical value'
    );
    assert(
      Math.abs(fx.integrate(-2.0, 2.0) / intfx) < 1.001,
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
  it('Ortho-Normalization', () => {
    let fx1 = new hamiltonian.BasisFunction('x + 2', 1.0, -1.0);
    let fx2 = new hamiltonian.BasisFunction('sin(x) + 2', 1.0, -1.0);
    let fx3 = new hamiltonian.BasisFunction('cos(x)', 1.0, -1.0);
    let basis = new hamiltonian.Basis([fx1, fx2, fx3]);
    basis.orthoNormalize();
    assert(basis.overlap(0, 1) <= 0.000001, 'Overlap should be 0 {01}');
    assert(basis.overlap(0, 2) <= 0.000001, 'Overlap should be 0 {02}');
    assert(basis.overlap(1, 2) <= 0.000001, 'Overlap should be 0 {12}');
    assert(Math.abs(basis.overlap(2, 2) - 1.0) <= 0.000001, 'Function not normalized');
    assert(Math.abs(basis.overlap(1, 1) - 1.0) <= 0.000001, 'Function not normalized');
    assert(Math.abs(basis.overlap(0, 0) - 1.0) <= 0.000001, 'Function not normalized');
  });
  it('2nd order Derivative Coupling', () => {
    let fx1 = new hamiltonian.BasisFunction('x + 2');
    let fx2 = new hamiltonian.BasisFunction('sin(x) + 2');
    let basis = new hamiltonian.Basis([fx1, fx2]);
    // eslint-disable-next-line new-cap
    assert(basis.Dij(1, 1) === -basis.set[1].kinetic(), 'Dij evaluation wrong');
  });
  it('Basis Normalization', () => {
    let fx1 = new hamiltonian.BasisFunction('x + 2');
    let fx2 = new hamiltonian.BasisFunction('sin(x) + 2');
    let basis = new hamiltonian.Basis([fx1, fx2]);
    basis.normalize();
    // eslint-disable-next-line new-cap
    assert(Math.abs(basis.overlap(0, 0) - 1.0) <= 0.001, 'Function not normalized');
    assert(Math.abs(basis.overlap(1, 1) - 1.0) <= 0.001, 'Function not normalized');
  });
});
describe('Hamiltonian Matrix', () => {
  it('Creation of a Hamiltonian Matrix', () => {
    let _;
    let fx1 = new hamiltonian.BasisFunction('sin(x)');
    let fx2 = new hamiltonian.BasisFunction('cos(x)');
    let basis = new hamiltonian.Basis([fx1, fx2]);
    basis.orthoNormalize();
    let H = hamiltonian.hamiltonian(basis, _, _, false);
    assert(H.length === basis.set.length, 'Hamiltonian Dimention not correct');
    let H2 = hamiltonian.hamiltonian(basis, _, _, true);
    assert(
      Math.abs(H[1][1] - H2[1][1]) <= 0.00000001,
      'Hamiltonian in orthogonal basis not correct'
    );
  });
});

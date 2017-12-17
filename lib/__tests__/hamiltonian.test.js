const assert = require('assert');
const hamiltonian = require('../index.js');

describe('hamiltonian', () => {
  it('has a test', () => {
    console.log(new hamiltonian.BasisFunction('x + 2'));
    assert(false, 'hamiltonian should have a test');
  });
});

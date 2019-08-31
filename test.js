
import {test} from 'tape';
import {orient2d} from './predicates.js';

test('orient2d', (t) => {
    const d = orient2d(0.9500000000002207, 0.9500000000002269, 18, 18, 16.8, 16.8);
    t.equal(d, -7.460698725481048e-15);
    t.end();
});

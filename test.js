
import {test} from 'tape';
import {orient2d} from './index.js';
import orient2dOld from 'robust-orientation';

test('orient2d', (t) => {
    t.ok(orient2d(0, 0, 1, 1, 0, 1) > 0, 'basic clockwise');
    t.ok(orient2d(0, 0, 0, 1, 1, 1) < 0, 'basic counter-clockwise');
    t.ok(orient2d(0, 0, 0.5, 0.5, 1, 1) === 0, 'basic collinear');

    const r = 0.95;
    const q = 18;
    const p = 16.8;
    const w = Math.pow(2, -43);

    for (let i = 0; i < 512; i++) {
        for (let j = 0; j < 512; j++) {
            const x = r + w * i / 512;
            const y = r + w * j / 512;

            const o = orient2d(x, y, q, q, p, p);
            const o2 = orient2dOld[3]([x, y], [p, p], [q, q]);

            if (Math.sign(o) !== Math.sign(o2)) {
                t.fail(`${x},${y}, ${q},${q}, ${p},${p}: ${o} vs ${o2}`);
            }
        }
    }
    t.pass('512x512 near-collinear');

    t.end();
});

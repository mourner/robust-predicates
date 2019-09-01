
import {test} from 'tape';
import orient2dOld from 'robust-orientation';
import nextafter from 'nextafter';

import {orient2d, orient3d, incircle} from './index.js';

test('orient2d', (t) => {
    t.ok(orient2d(0, 0, 1, 1, 0, 1) > 0, 'clockwise');
    t.ok(orient2d(0, 0, 0, 1, 1, 1) < 0, 'counter-clockwise');
    t.ok(orient2d(0, 0, 0.5, 0.5, 1, 1) === 0, 'collinear');

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

test('incircle', (t) => {
    t.ok(incircle(0, -1, 1, 0, 0, 1, -0.5, 0) > 0, 'inside');
    t.ok(incircle(0, -1, 1, 0, 0, 1, -1, 0) === 0, 'on circle');
    t.ok(incircle(0, -1, 1, 0, 0, 1, -1.5, 0) < 0, 'outside');

    const a = nextafter(-1, 0);
    const b = nextafter(-1, -2);

    t.ok(incircle(1, 0, 0, 1, -1, 0, 0, a) > 0, 'near inside');
    t.ok(incircle(1, 0, 0, 1, -1, 0, 0, b) < 0, 'near outside');

    t.end();
});

test('orient3d', (t) => {
    t.ok(orient3d(
        0, 0, 0,
        0, 1, 0,
        1, 0, 0,
        0, 0, 1
    ) > 0, 'below');

    t.ok(orient3d(
        0, 0, 0,
        0, 1, 0,
        1, 0, 0,
        0, 0, -1
    ) < 0, 'above');

    t.ok(orient3d(
        0, 0, 0,
        0, 1, 0,
        1, 0, 0,
        0, 0, 0
    ) === 0, 'coplanar');

    const a = nextafter(0, 1);
    const b = nextafter(0, -1);

    t.ok(orient3d(
        0, 0, 0,
        0, 1, 0,
        1, 0, 0,
        0, 0, a
    ) > 0, 'near below');

    t.ok(orient3d(
        0, 0, 0,
        0, 1, 0,
        1, 0, 0,
        0, 0, b
    ) < 0, 'near above');

    t.end();
});

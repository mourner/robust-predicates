
import {test} from 'tape';
import robustOrientation from 'robust-orientation';
import nextafter from 'nextafter';

import {
    orient2d, orient2dfast,
    orient3d, orient3dfast,
    incircle, incirclefast,
    insphere, inspherefast
} from './index.js';

test('orient2d', (t) => {
    t.ok(orient2d(0, 0, 1, 1, 0, 1) < 0, 'clockwise');
    t.ok(orient2d(0, 0, 0, 1, 1, 1) > 0, 'counterclockwise');
    t.ok(orient2d(0, 0, 0.5, 0.5, 1, 1) === 0, 'collinear');

    const r = 0.95;
    const q = 18;
    const p = 16.8;
    const w = Math.pow(2, -43);

    for (let i = 0; i < 128; i++) {
        for (let j = 0; j < 128; j++) {
            const x = r + w * i / 128;
            const y = r + w * j / 128;

            const o = orient2d(x, y, q, q, p, p);
            const o2 = robustOrientation[3]([x, y], [q, q], [p, p]);

            if (Math.sign(o) !== Math.sign(o2)) {
                t.fail(`${x},${y}, ${q},${q}, ${p},${p}: ${o} vs ${o2}`);
            }
        }
    }
    t.pass('512x512 near-collinear');

    t.end();
});

test('orient2dfast', (t) => {
    t.ok(orient2dfast(0, 0, 1, 1, 0, 1) < 0, 'counterclockwise');
    t.ok(orient2dfast(0, 0, 0, 1, 1, 1) > 0, 'clockwise');
    t.ok(orient2dfast(0, 0, 0.5, 0.5, 1, 1) === 0, 'collinear');
    t.end();
});

test('incircle', (t) => {
    t.ok(incircle(0, -1, 0, 1, 1, 0, -0.5, 0) < 0, 'inside');
    t.ok(incircle(0, -1, 1, 0, 0, 1, -1, 0) === 0, 'on circle');
    t.ok(incircle(0, -1, 0, 1, 1, 0, -1.5, 0) > 0, 'outside');

    const a = nextafter(-1, 0);
    const b = nextafter(-1, -2);

    t.ok(incircle(1, 0, -1, 0, 0, 1, 0, a) < 0, 'near inside');
    t.ok(incircle(1, 0, -1, 0, 0, 1, 0, b) > 0, 'near outside');

    let x = 1e-64;
    for (let i = 0; i < 128; i++) {
        if (incircle(0, x, -x, -x, x, -x, 0, 0) <= 0) t.fail(`incircle test ${x}, outside`);
        if (incircle(0, x, -x, -x, x, -x, 0, 2 * x) >= 0) t.fail(`incircle test ${x}, inside`);
        if (incircle(0, x, -x, -x, x, -x, 0, x) !== 0) t.fail(`incircle test ${x}, cocircular`);
        x *= 10;
    }
    t.pass(`${128 * 3} incircle tests`);

    t.end();
});

test('incirclefast', (t) => {
    t.ok(incirclefast(0, -1, 0, 1, 1, 0, -0.5, 0) < 0, 'inside');
    t.ok(incirclefast(0, -1, 0, 1, 1, 0, -1, 0) === 0, 'on circle');
    t.ok(incirclefast(0, -1, 0, 1, 1, 0, -1.5, 0) > 0, 'outside');
    t.end();
});

test('orient3d', (t) => {
    t.ok(orient3d(
        0, 0, 0,
        0, 1, 0,
        1, 0, 0,
        0, 0, 1
    ) > 0, 'above');

    t.ok(orient3d(
        0, 0, 0,
        0, 1, 0,
        1, 0, 0,
        0, 0, -1
    ) < 0, 'below');

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
    ) > 0, 'near above');

    t.ok(orient3d(
        0, 0, 0,
        0, 1, 0,
        1, 0, 0,
        0, 0, b
    ) < 0, 'near below');

    t.end();
});

test('orient3dfast', (t) => {
    t.ok(orient3dfast(
        0, 0, 0,
        0, 1, 0,
        1, 0, 0,
        0, 0, 1
    ) > 0, 'above');

    t.ok(orient3dfast(
        0, 0, 0,
        0, 1, 0,
        1, 0, 0,
        0, 0, -1
    ) < 0, 'below');

    t.ok(orient3dfast(
        0, 0, 0,
        0, 1, 0,
        1, 0, 0,
        0, 0, 0
    ) === 0, 'coplanar');

    t.end();
});

test('insphere', (t) => {
    t.ok(insphere(
        1, 0, 0,
        0, -1, 0,
        0, 1, 0,
        0, 0, 1,
        0, 0, 0
    ) < 0, 'inside');

    t.ok(insphere(
        1, 0, 0,
        0, -1, 0,
        0, 1, 0,
        0, 0, 1,
        0, 0, 2
    ) > 0, 'outside');

    t.ok(insphere(
        1, 0, 0,
        0, -1, 0,
        0, 1, 0,
        0, 0, 1,
        0, 0, -1
    ) === 0, 'cospherical');

    const a = nextafter(-1, 0);
    const b = nextafter(-1, -2);

    t.ok(insphere(
        1, 0, 0,
        0, -1, 0,
        0, 1, 0,
        0, 0, 1,
        0, 0, a
    ) < 0, 'near inside');

    t.ok(insphere(
        1, 0, 0,
        0, -1, 0,
        0, 1, 0,
        0, 0, 1,
        0, 0, b
    ) > 0, 'near outside');

    t.end();
});

test('inspherefast', (t) => {
    t.ok(inspherefast(
        1, 0, 0,
        0, -1, 0,
        0, 1, 0,
        0, 0, 1,
        0, 0, 0
    ) < 0, 'inside');

    t.ok(inspherefast(
        1, 0, 0,
        0, -1, 0,
        0, 1, 0,
        0, 0, 1,
        0, 0, 2
    ) > 0, 'outside');

    t.ok(inspherefast(
        1, 0, 0,
        0, -1, 0,
        0, 1, 0,
        0, 0, 1,
        0, 0, -1
    ) === 0, 'cospherical');

    t.end();
});

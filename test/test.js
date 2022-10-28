
import fs from 'fs';
import test from 'node:test';
import assert from 'node:assert/strict';
import robustOrientation from 'robust-orientation';
import nextafter from 'nextafter';

import {
    orient2d, orient2dfast,
    orient3d, orient3dfast,
    incircle, incirclefast,
    insphere, inspherefast
} from '../index.js';

test('orient2d', () => {
    assert.ok(orient2d(0, 0, 1, 1, 0, 1) < 0, 'clockwise');
    assert.ok(orient2d(0, 0, 0, 1, 1, 1) > 0, 'counterclockwise');
    assert.ok(orient2d(0, 0, 0.5, 0.5, 1, 1) === 0, 'collinear');

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
                assert.fail(`${x},${y}, ${q},${q}, ${p},${p}: ${o} vs ${o2}`);
            }
        }
    }
    // 512x512 near-collinear

    const lines = fs.readFileSync(new URL('./fixtures/orient2d.txt', import.meta.url), 'utf8').trim().split(/\r?\n/);
    for (const line of lines) {
        const [, ax, ay, bx, by, cx, cy, sign] = line.split(' ').map(Number);
        const result = orient2d(ax, ay, bx, by, cx, cy);
        if (Math.sign(result) !== -sign) {
            assert.fail(`${line}: ${result} vs ${-sign}`);
        }
    }
    // 1000 hard fixtures
});

test('orient2dfast', () => {
    assert.ok(orient2dfast(0, 0, 1, 1, 0, 1) < 0, 'counterclockwise');
    assert.ok(orient2dfast(0, 0, 0, 1, 1, 1) > 0, 'clockwise');
    assert.ok(orient2dfast(0, 0, 0.5, 0.5, 1, 1) === 0, 'collinear');
});

test('incircle', () => {
    assert.ok(incircle(0, -1, 0, 1, 1, 0, -0.5, 0) < 0, 'inside');
    assert.ok(incircle(0, -1, 1, 0, 0, 1, -1, 0) === 0, 'on circle');
    assert.ok(incircle(0, -1, 0, 1, 1, 0, -1.5, 0) > 0, 'outside');

    const a = nextafter(-1, 0);
    const b = nextafter(-1, -2);

    assert.ok(incircle(1, 0, -1, 0, 0, 1, 0, a) < 0, 'near inside');
    assert.ok(incircle(1, 0, -1, 0, 0, 1, 0, b) > 0, 'near outside');

    let x = 1e-64;
    for (let i = 0; i < 128; i++) {
        if (incircle(0, x, -x, -x, x, -x, 0, 0) <= 0) assert.fail(`incircle test ${x}, outside`);
        if (incircle(0, x, -x, -x, x, -x, 0, 2 * x) >= 0) assert.fail(`incircle test ${x}, inside`);
        if (incircle(0, x, -x, -x, x, -x, 0, x) !== 0) assert.fail(`incircle test ${x}, cocircular`);
        x *= 10;
    }
    // 384 incircle tests

    const lines = fs.readFileSync(new URL('./fixtures/incircle.txt', import.meta.url), 'utf8').trim().split(/\r?\n/);
    for (const line of lines) {
        const [, ax, ay, bx, by, cx, cy, dx, dy, sign] = line.split(' ').map(Number);
        const result = incircle(ax, ay, bx, by, cx, cy, dx, dy);
        if (Math.sign(result) !== sign) {
            assert.fail(`${line}: ${result} vs ${sign}`);
        }
    }
    // 1000 hard fixtures
});

test('incirclefast', () => {
    assert.ok(incirclefast(0, -1, 0, 1, 1, 0, -0.5, 0) < 0, 'inside');
    assert.ok(incirclefast(0, -1, 0, 1, 1, 0, -1, 0) === 0, 'on circle');
    assert.ok(incirclefast(0, -1, 0, 1, 1, 0, -1.5, 0) > 0, 'outside');
});

test('orient3d', () => {
    assert.ok(orient3d(
        0, 0, 0,
        0, 1, 0,
        1, 0, 0,
        0, 0, 1
    ) > 0, 'above');

    assert.ok(orient3d(
        0, 0, 0,
        0, 1, 0,
        1, 0, 0,
        0, 0, -1
    ) < 0, 'below');

    assert.ok(orient3d(
        0, 0, 0,
        0, 1, 0,
        1, 0, 0,
        0, 0, 0
    ) === 0, 'coplanar');

    const a = nextafter(0, 1);
    const b = nextafter(0, -1);

    assert.ok(orient3d(
        0, 0, 0,
        0, 1, 0,
        1, 0, 0,
        0, 0, a
    ) > 0, 'near above');

    assert.ok(orient3d(
        0, 0, 0,
        0, 1, 0,
        1, 0, 0,
        0, 0, b
    ) < 0, 'near below');

    const lines = fs.readFileSync(new URL('./fixtures/orient3d.txt', import.meta.url), 'utf8').trim().split(/\r?\n/);
    for (const line of lines) {
        const [, ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz, sign] = line.split(' ').map(Number);
        const result = orient3d(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz);
        if (Math.sign(result) !== sign) assert.fail(`${line}: ${result} vs ${sign}`);
        if (Math.sign(result) !== Math.sign(orient3d(dx, dy, dz, bx, by, bz, ax, ay, az, cx, cy, cz))) assert.fail('symmetry');
    }
    // 1000 hard fixtures

    const tol = 5.0e-14;

    for (let i = 0; i < 1000; i++) {
        const ax = 0.5 + tol * Math.random();
        const ay = 0.5 + tol * Math.random();
        const az = 0.5 + tol * Math.random();
        const b = 12, c = 24, d = 48;
        if (orient3d(b, b, b, c, c, c, d, d, d, ax, ay, az) !== 0) assert.fail('degenerate');
        if (orient3d(c, c, c, d, d, d, ax, ay, az, b, b, b) !== 0) assert.fail('degenerate');
    }
    // 1000 degenerate cases
});

test('orient3dfast', () => {
    assert.ok(orient3dfast(
        0, 0, 0,
        0, 1, 0,
        1, 0, 0,
        0, 0, 1
    ) > 0, 'above');

    assert.ok(orient3dfast(
        0, 0, 0,
        0, 1, 0,
        1, 0, 0,
        0, 0, -1
    ) < 0, 'below');

    assert.ok(orient3dfast(
        0, 0, 0,
        0, 1, 0,
        1, 0, 0,
        0, 0, 0
    ) === 0, 'coplanar');
});

test('insphere', () => {
    assert.ok(insphere(
        1, 0, 0,
        0, -1, 0,
        0, 1, 0,
        0, 0, 1,
        0, 0, 0
    ) < 0, 'inside');

    assert.ok(insphere(
        1, 0, 0,
        0, -1, 0,
        0, 1, 0,
        0, 0, 1,
        0, 0, 2
    ) > 0, 'outside');

    assert.ok(insphere(
        1, 0, 0,
        0, -1, 0,
        0, 1, 0,
        0, 0, 1,
        0, 0, -1
    ) === 0, 'cospherical');

    const a = nextafter(-1, 0);
    const b = nextafter(-1, -2);

    assert.ok(insphere(
        1, 0, 0,
        0, -1, 0,
        0, 1, 0,
        0, 0, 1,
        0, 0, a
    ) < 0, 'near inside');

    assert.ok(insphere(
        1, 0, 0,
        0, -1, 0,
        0, 1, 0,
        0, 0, 1,
        0, 0, b
    ) > 0, 'near outside');

    const lines = fs.readFileSync(new URL('./fixtures/insphere.txt', import.meta.url), 'utf8').trim().split(/\r?\n/);
    for (const line of lines) {
        const [, ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz, ex, ey, ez, sign] = line.split(' ').map(Number);
        const result = insphere(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz, ex, ey, ez);
        if (Math.sign(result) !== -sign) {
            assert.fail(`${line}: ${result} vs ${-sign}`);
        }
    }
    // 1000 hard fixtures
});

test('inspherefast', () => {
    assert.ok(inspherefast(
        1, 0, 0,
        0, -1, 0,
        0, 1, 0,
        0, 0, 1,
        0, 0, 0
    ) < 0, 'inside');

    assert.ok(inspherefast(
        1, 0, 0,
        0, -1, 0,
        0, 1, 0,
        0, 0, 1,
        0, 0, 2
    ) > 0, 'outside');

    assert.ok(inspherefast(
        1, 0, 0,
        0, -1, 0,
        0, 1, 0,
        0, 0, 1,
        0, 0, -1
    ) === 0, 'cospherical');
});

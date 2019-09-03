
const fs = require('fs');
const path = require('path');
const robustOrientation = require('robust-orientation');
const robustInSphere = require('robust-in-sphere');

const {orient2d, orient3d, incircle, insphere} = require('./umd/predicates.js');

{
    const lines = fs.readFileSync(path.join(__dirname, 'test/fixtures/orient2d.txt'), 'utf8').trim().split(/\r?\n/);
    const coords = new Float64Array(lines.length * 6);
    const points = [];
    let i = 0;
    for (const line of lines) {
        const [, ax, ay, bx, by, cx, cy] = line.split(' ').map(Number);
        coords[i++] = ax;
        coords[i++] = ay;
        coords[i++] = bx;
        coords[i++] = by;
        coords[i++] = cx;
        coords[i++] = cy;
        points.push([ax, ay], [bx, by], [cx, cy]);
    }
    const N = 100000;
    const id = `${N * lines.length} x orient2d `;
    console.time(`${id}(robust-predicates)`);
    for (let k = 0; k < N; k++) {
        for (let i = 0; i < coords.length; i += 6) orient2d(
            coords[i + 0], coords[i + 1],
            coords[i + 2], coords[i + 3],
            coords[i + 4], coords[i + 5]
        );
    }
    console.timeEnd(`${id}(robust-predicates)`);

    console.time(`${id}(robust-orientation)`);
    for (let k = 0; k < N; k++) {
        for (let i = 0; i < points.length; i += 3) robustOrientation[3](
            points[i + 0],
            points[i + 1],
            points[i + 2]
        );
    }
    console.timeEnd(`${id}(robust-orientation)`);
}

{
    const lines = fs.readFileSync(path.join(__dirname, 'test/fixtures/orient3d.txt'), 'utf8').trim().split(/\r?\n/);
    const coords = new Float64Array(lines.length * 12);
    const points = [];
    let i = 0;
    for (const line of lines) {
        const [, ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz] = line.split(' ').map(Number);
        coords[i++] = ax;
        coords[i++] = ay;
        coords[i++] = az;
        coords[i++] = bx;
        coords[i++] = by;
        coords[i++] = bz;
        coords[i++] = cx;
        coords[i++] = cy;
        coords[i++] = cz;
        coords[i++] = dx;
        coords[i++] = dy;
        coords[i++] = dz;
        points.push([ax, ay, az], [bx, by, bz], [cx, cy, cz], [dx, dy, dz]);
    }
    const N = 10000;
    const id = `${N * lines.length} x orient3d `;
    console.time(`${id}(robust-predicates)`);
    for (let k = 0; k < N; k++) {
        for (let i = 0; i < coords.length; i += 12) orient3d(
            coords[i + 0], coords[i + 1], coords[i + 2],
            coords[i + 3], coords[i + 4], coords[i + 5],
            coords[i + 6], coords[i + 7], coords[i + 8],
            coords[i + 9], coords[i + 10], coords[i + 11]
        );
    }
    console.timeEnd(`${id}(robust-predicates)`);

    console.time(`${id}(robust-orientation)`);
    for (let k = 0; k < N; k++) {
        for (let i = 0; i < points.length; i += 4) robustOrientation[4](
            points[i + 0],
            points[i + 1],
            points[i + 2],
            points[i + 3]
        );
    }
    console.timeEnd(`${id}(robust-orientation)`);
}

{
    const lines = fs.readFileSync(path.join(__dirname, 'test/fixtures/incircle.txt'), 'utf8').trim().split(/\r?\n/);
    const coords = new Float64Array(lines.length * 8);
    const points = [];
    let i = 0;
    for (const line of lines) {
        const [, ax, ay, bx, by, cx, cy, dx, dy] = line.split(' ').map(Number);
        coords[i++] = ax;
        coords[i++] = ay;
        coords[i++] = bx;
        coords[i++] = by;
        coords[i++] = cx;
        coords[i++] = cy;
        coords[i++] = dx;
        coords[i++] = dy;
        points.push([ax, ay], [bx, by], [cx, cy], [dx, dy]);
    }
    const N = 100;
    const id = `${N * lines.length} x incircle `;
    console.time(`${id}(robust-predicates)`);
    for (let k = 0; k < N; k++) {
        for (let i = 0; i < coords.length; i += 8) incircle(
            coords[i + 0], coords[i + 1],
            coords[i + 2], coords[i + 3],
            coords[i + 4], coords[i + 5],
            coords[i + 6], coords[i + 7]
        );
    }
    console.timeEnd(`${id}(robust-predicates)`);

    console.time(`${id}(robust-in-sphere)`);
    for (let k = 0; k < N; k++) {
        for (let i = 0; i < points.length; i += 4) robustInSphere[4](
            points[i + 0],
            points[i + 1],
            points[i + 2],
            points[i + 3]
        );
    }
    console.timeEnd(`${id}(robust-in-sphere)`);
}

{
    const lines = fs.readFileSync(path.join(__dirname, 'test/fixtures/insphere.txt'), 'utf8').trim().split(/\r?\n/);
    const coords = new Float64Array(lines.length * 15);
    const points = [];
    let i = 0;
    for (const line of lines) {
        const [, ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz, ex, ey, ez] = line.split(' ').map(Number);
        coords[i++] = ax;
        coords[i++] = ay;
        coords[i++] = az;
        coords[i++] = bx;
        coords[i++] = by;
        coords[i++] = bz;
        coords[i++] = cx;
        coords[i++] = cy;
        coords[i++] = cz;
        coords[i++] = dx;
        coords[i++] = dy;
        coords[i++] = dz;
        coords[i++] = ex;
        coords[i++] = ey;
        coords[i++] = ez;
        points.push([ax, ay, az], [bx, by, bz], [cx, cy, cz], [dx, dy, dz], [ex, ey, ez]);
    }
    const N = 10;
    const id = `${N * lines.length} x insphere `;
    console.time(`${id}(robust-predicates)`);
    for (let k = 0; k < N; k++) {
        for (let i = 0; i < coords.length; i += 15) insphere(
            coords[i + 0], coords[i + 1], coords[i + 2],
            coords[i + 3], coords[i + 4], coords[i + 5],
            coords[i + 6], coords[i + 7], coords[i + 8],
            coords[i + 9], coords[i + 10], coords[i + 11],
            coords[i + 12], coords[i + 13], coords[i + 14]
        );
    }
    console.timeEnd(`${id}(robust-predicates)`);

    console.time(`${id}(robust-in-sphere)`);
    for (let k = 0; k < N; k++) {
        for (let i = 0; i < points.length; i += 5) robustInSphere[5](
            points[i + 0],
            points[i + 1],
            points[i + 2],
            points[i + 3],
            points[i + 4]
        );
    }
    console.timeEnd(`${id}(robust-in-sphere)`);
}

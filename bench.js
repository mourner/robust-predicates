
import {orient2d} from './predicates.js';
import orientOld from 'robust-orientation';

const r = 0.95;
const q = [18, 18];
const p = [16.8, 16.8];
const w = Math.pow(2, -43);
const points = [];
const size = 1024;

for (let i = 0; i < size; i++) {
    for (let j = 0; j < size; j++) {
        const x = r + w * i / size;
        const y = r + w * j / size;
        points.push([x, y]);
    }
}

console.time('robust-predicates orient2d');
for (let r of points) orient2d(r[0], r[1], q[0], q[1], p[0], p[1]);
console.timeEnd('robust-predicates orient2d');

console.time('robust-orientation');
for (let r of points) orientOld[3](r, p, q);
console.timeEnd('robust-orientation');

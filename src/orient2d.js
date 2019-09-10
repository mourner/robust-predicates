import {epsilon, splitter, resulterrbound, estimate, vec, sum} from './util.js';

const ccwerrboundA = (3 + 16 * epsilon) * epsilon;
const ccwerrboundB = (2 + 12 * epsilon) * epsilon;
const ccwerrboundC = (9 + 64 * epsilon) * epsilon * epsilon;

const B = vec(4);
const C1 = vec(8);
const C2 = vec(12);
const D = vec(16);
const u = vec(4);

function orient2dadapt(ax, ay, bx, by, cx, cy, detsum) {
    let acxtail, acytail, bcxtail, bcytail;
    let bvirt, c, ahi, alo, bhi, blo, _i, _j, _0, s1, s0, t1, t0, u3;

    const acx = ax - cx;
    const bcx = bx - cx;
    const acy = ay - cy;
    const bcy = by - cy;

    $Cross_Product(acx, bcx, acy, bcy, B);

    let det = estimate(4, B);
    let errbound = ccwerrboundB * detsum;
    if (det >= errbound || -det >= errbound) {
        return det;
    }

    $Two_Diff_Tail(ax, cx, acx, acxtail);
    $Two_Diff_Tail(bx, cx, bcx, bcxtail);
    $Two_Diff_Tail(ay, cy, acy, acytail);
    $Two_Diff_Tail(by, cy, bcy, bcytail);

    if (acxtail === 0 && acytail === 0 && bcxtail === 0 && bcytail === 0) {
        return det;
    }

    errbound = ccwerrboundC * detsum + resulterrbound * Math.abs(det);
    det += (acx * bcytail + bcy * acxtail) - (acy * bcxtail + bcx * acytail);
    if (det >= errbound || -det >= errbound) return det;

    $Cross_Product(acxtail, bcx, acytail, bcy, u);
    const C1len = sum(4, B, 4, u, C1);

    $Cross_Product(acx, bcxtail, acy, bcytail, u);
    const C2len = sum(C1len, C1, 4, u, C2);

    $Cross_Product(acxtail, bcxtail, acytail, bcytail, u);
    const Dlen = sum(C2len, C2, 4, u, D);

    return D[Dlen - 1];
}

export function orient2d(ax, ay, bx, by, cx, cy) {
    const detleft = (ay - cy) * (bx - cx);
    const detright = (ax - cx) * (by - cy);
    const det = detleft - detright;

    if (detleft === 0 || detright === 0 || (detleft > 0) !== (detright > 0)) return det;

    const detsum = Math.abs(detleft + detright);
    if (Math.abs(det) >= ccwerrboundA * detsum) return det;

    return -orient2dadapt(ax, ay, bx, by, cx, cy, detsum);
}

export function orient2dfast(ax, ay, bx, by, cx, cy) {
    return (ay - cy) * (bx - cx) - (ax - cx) * (by - cy);
}

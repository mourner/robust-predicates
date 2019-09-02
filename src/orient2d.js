import {
    epsilon, splitter, resulterrbound, estimate, vec,
    fast_expansion_sum_zeroelim
} from './util.js';

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
    let detleft, detright, detlefttail, detrighttail;
    let B3, u3, s1, t1, s0, t0;
    let bvirt, c, ahi, alo, bhi, blo, _i, _j, _0;

    const acx = ax - cx;
    const bcx = bx - cx;
    const acy = ay - cy;
    const bcy = by - cy;

    $Two_Product(acx, bcy, detleft, detlefttail);
    $Two_Product(acy, bcx, detright, detrighttail);

    $Two_Two_Diff(detleft, detlefttail, detright, detrighttail, B3, B[2], B[1], B[0]);
    B[3] = B3;

    let det = estimate(4, B);
    let errbound = ccwerrboundB * detsum;
    if (det >= errbound || -det >= errbound) return det;

    $Two_Diff_Tail(ax, cx, acx, acxtail);
    $Two_Diff_Tail(bx, cx, bcx, bcxtail);
    $Two_Diff_Tail(ay, cy, acy, acytail);
    $Two_Diff_Tail(by, cy, bcy, bcytail);

    if (acxtail === 0 && acytail === 0 && bcxtail === 0 && bcytail === 0) return det;

    errbound = ccwerrboundC * detsum + resulterrbound * Math.abs(det);
    det += (acx * bcytail + bcy * acxtail) - (acy * bcxtail + bcx * acytail);
    if (det >= errbound || -det >= errbound) {
        return det;
    }

    $Two_Product(acxtail, bcy, s1, s0);
    $Two_Product(acytail, bcx, t1, t0);
    $Two_Two_Diff(s1, s0, t1, t0, u3, u[2], u[1], u[0]);
    u[3] = u3;
    const C1length = fast_expansion_sum_zeroelim(4, B, 4, u, C1);

    $Two_Product(acx, bcytail, s1, s0);
    $Two_Product(acy, bcxtail, t1, t0);
    $Two_Two_Diff(s1, s0, t1, t0, u3, u[2], u[1], u[0]);
    u[3] = u3;
    const C2length = fast_expansion_sum_zeroelim(C1length, C1, 4, u, C2);

    $Two_Product(acxtail, bcytail, s1, s0);
    $Two_Product(acytail, bcxtail, t1, t0);
    $Two_Two_Diff(s1, s0, t1, t0, u3, u[2], u[1], u[0]);
    u[3] = u3;
    const Dlength = fast_expansion_sum_zeroelim(C2length, C2, 4, u, D);

    return D[Dlength - 1];
}

export function orient2d(ax, ay, bx, by, cx, cy) {
    const detleft = (ay - cy) * (bx - cx);
    const detright = (ax - cx) * (by - cy);
    const det = detleft - detright;
    let detsum;

    if (detleft > 0) {
        if (detright <= 0) {
            return det;
        }
        detsum = detleft + detright;

    } else if (detleft < 0) {
        if (detright >= 0) {
            return det;
        }
        detsum = -detleft - detright;

    } else {
        return det;
    }

    const errbound = ccwerrboundA * detsum;
    if (det >= errbound || -det >= errbound) {
        return det;
    }

    return -orient2dadapt(ax, ay, bx, by, cx, cy, detsum);
}

export function orient2dfast(ax, ay, bx, by, cx, cy) {
    return (ay - cy) * (bx - cx) - (ax - cx) * (by - cy);
}

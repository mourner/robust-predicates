import {epsilon, splitter, resulterrbound, estimate, vec, sum, scale} from './util.js';

const o3derrboundA = (7 + 56 * epsilon) * epsilon;
const o3derrboundB = (3 + 28 * epsilon) * epsilon;
const o3derrboundC = (26 + 288 * epsilon) * epsilon * epsilon;

const bc = vec(4);
const ca = vec(4);
const ab = vec(4);
const at_b = vec(4);
const at_c = vec(4);
const bt_c = vec(4);
const bt_a = vec(4);
const ct_a = vec(4);
const ct_b = vec(4);
const bct = vec(8);
const cat = vec(8);
const abt = vec(8);
const u = vec(4);

const _8 = vec(8);
const _8b = vec(8);
const _16 = vec(8);
const _12 = vec(12);

let fin = vec(192);
let fin2 = vec(192);

function finadd(finlen, alen, a) {
    finlen = sum(finlen, fin, alen, a, fin2);
    const tmp = fin; fin = fin2; fin2 = tmp;
    return finlen;
}

function tailinit(xtail, ytail, ax, ay, bx, by, a, b) {
    let bvirt, c, ahi, alo, bhi, blo, _i, _j, _k, _0, s1, s0, t1, t0, u3, negate;
    if (xtail === 0) {
        if (ytail === 0) {
            a[0] = 0;
            b[0] = 0;
            return 1;
        }
        negate = -ytail;
        $Two_Product(negate, ax, s1, a[0]);
        a[1] = s1;
        $Two_Product(ytail, bx, s1, b[0]);
        b[1] = s1;
        return 2;
    }
    if (ytail === 0) {
        $Two_Product(xtail, ay, s1, a[0]);
        a[1] = s1;
        negate = -xtail;
        $Two_Product(negate, by, s1, b[0]);
        b[1] = s1;
        return 2;
    }
    $Cross_Product(xtail, ax, ytail, ay, a);
    $Cross_Product(ytail, by, xtail, bx, b);
    return 4;
}

function tailadd(finlen, a, b, k, z) {
    let bvirt, c, ahi, alo, bhi, blo, _i, _j, _k, _0, s1, s0, u3;
    $Two_Product(a, b, s1, s0);
    $Two_One_Product(s1, s0, k, u);
    finlen = finadd(finlen, 4, u);
    if (z !== 0) {
        $Two_One_Product(s1, s0, z, u);
        finlen = finadd(finlen, 4, u);
    }
    return finlen;
}

function orient3dadapt(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz, permanent) {
    let finlen;
    let adxtail, bdxtail, cdxtail;
    let adytail, bdytail, cdytail;
    let adztail, bdztail, cdztail;
    let bvirt, c, ahi, alo, bhi, blo, _i, _j, _k, _0, s1, s0, t1, t0, u3;

    const adx = ax - dx;
    const bdx = bx - dx;
    const cdx = cx - dx;
    const ady = ay - dy;
    const bdy = by - dy;
    const cdy = cy - dy;
    const adz = az - dz;
    const bdz = bz - dz;
    const cdz = cz - dz;

    $Cross_Product(bdx, bdy, cdx, cdy, bc);
    $Cross_Product(cdx, cdy, adx, ady, ca);
    $Cross_Product(adx, ady, bdx, bdy, ab);

    finlen = sum(
        sum(
            scale(4, bc, adz, _8), _8,
            scale(4, ca, bdz, _8b), _8b, _16), _16,
        scale(4, ab, cdz, _8), _8, fin);

    let det = estimate(finlen, fin);
    let errbound = o3derrboundB * permanent;
    if (det >= errbound || -det >= errbound) {
        return det;
    }

    $Two_Diff_Tail(ax, dx, adx, adxtail);
    $Two_Diff_Tail(bx, dx, bdx, bdxtail);
    $Two_Diff_Tail(cx, dx, cdx, cdxtail);
    $Two_Diff_Tail(ay, dy, ady, adytail);
    $Two_Diff_Tail(by, dy, bdy, bdytail);
    $Two_Diff_Tail(cy, dy, cdy, cdytail);
    $Two_Diff_Tail(az, dz, adz, adztail);
    $Two_Diff_Tail(bz, dz, bdz, bdztail);
    $Two_Diff_Tail(cz, dz, cdz, cdztail);

    if (adxtail === 0 && bdxtail === 0 && cdxtail === 0 &&
        adytail === 0 && bdytail === 0 && cdytail === 0 &&
        adztail === 0 && bdztail === 0 && cdztail === 0) {
        return det;
    }

    errbound = o3derrboundC * permanent + resulterrbound * Math.abs(det);
    det +=
        adz * (bdx * cdytail + cdy * bdxtail - (bdy * cdxtail + cdx * bdytail)) + adztail * (bdx * cdy - bdy * cdx) +
        bdz * (cdx * adytail + ady * cdxtail - (cdy * adxtail + adx * cdytail)) + bdztail * (cdx * ady - cdy * adx) +
        cdz * (adx * bdytail + bdy * adxtail - (ady * bdxtail + bdx * adytail)) + cdztail * (adx * bdy - ady * bdx);
    if (det >= errbound || -det >= errbound) {
        return det;
    }

    const at_len = tailinit(adxtail, adytail, bdx, bdy, cdx, cdy, at_b, at_c);
    const bt_len = tailinit(bdxtail, bdytail, cdx, cdy, adx, ady, bt_c, bt_a);
    const ct_len = tailinit(cdxtail, cdytail, adx, ady, bdx, bdy, ct_a, ct_b);

    const bctlen = sum(bt_len, bt_c, ct_len, ct_b, bct);
    finlen = finadd(finlen, scale(bctlen, bct, adz, _16), _16);

    const catlen = sum(ct_len, ct_a, at_len, at_c, cat);
    finlen = finadd(finlen, scale(catlen, cat, bdz, _16), _16);

    const abtlen = sum(at_len, at_b, bt_len, bt_a, abt);
    finlen = finadd(finlen, scale(abtlen, abt, cdz, _16), _16);

    if (adztail !== 0) {
        finlen = finadd(finlen, scale(4, bc, adztail, _12), _12);
        finlen = finadd(finlen, scale(bctlen, bct, adztail, _16), _16);
    }
    if (bdztail !== 0) {
        finlen = finadd(finlen, scale(4, ca, bdztail, _12), _12);
        finlen = finadd(finlen, scale(catlen, cat, bdztail, _16), _16);
    }
    if (cdztail !== 0) {
        finlen = finadd(finlen, scale(4, ab, cdztail, _12), _12);
        finlen = finadd(finlen, scale(abtlen, abt, cdztail, _16), _16);
    }

    if (adxtail !== 0) {
        if (bdytail !== 0) {
            finlen = tailadd(finlen, adxtail, bdytail, cdz, cdztail);
        }
        if (cdytail !== 0) {
            finlen = tailadd(finlen, -adxtail, cdytail, bdz, bdztail);
        }
    }
    if (bdxtail !== 0) {
        if (cdytail !== 0) {
            finlen = tailadd(finlen, bdxtail, cdytail, adz, adztail);
        }
        if (adytail !== 0) {
            finlen = tailadd(finlen, -bdxtail, adytail, cdz, cdztail);
        }
    }
    if (cdxtail !== 0) {
        if (adytail !== 0) {
            finlen = tailadd(finlen, cdxtail, adytail, bdz, bdztail);
        }
        if (bdytail !== 0) {
            finlen = tailadd(finlen, -cdxtail, bdytail, adz, adztail);
        }
    }

    return fin[finlen - 1];
}

export function orient3d(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz) {
    const adx = ax - dx;
    const bdx = bx - dx;
    const cdx = cx - dx;
    const ady = ay - dy;
    const bdy = by - dy;
    const cdy = cy - dy;
    const adz = az - dz;
    const bdz = bz - dz;
    const cdz = cz - dz;

    const bdxcdy = bdx * cdy;
    const cdxbdy = cdx * bdy;

    const cdxady = cdx * ady;
    const adxcdy = adx * cdy;

    const adxbdy = adx * bdy;
    const bdxady = bdx * ady;

    const det =
        adz * (bdxcdy - cdxbdy) +
        bdz * (cdxady - adxcdy) +
        cdz * (adxbdy - bdxady);

    const permanent =
        (Math.abs(bdxcdy) + Math.abs(cdxbdy)) * Math.abs(adz) +
        (Math.abs(cdxady) + Math.abs(adxcdy)) * Math.abs(bdz) +
        (Math.abs(adxbdy) + Math.abs(bdxady)) * Math.abs(cdz);

    const errbound = o3derrboundA * permanent;
    if (det > errbound || -det > errbound) {
        return det;
    }

    return orient3dadapt(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz, permanent);
}

export function orient3dfast(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz) {
    const adx = ax - dx;
    const bdx = bx - dx;
    const cdx = cx - dx;
    const ady = ay - dy;
    const bdy = by - dy;
    const cdy = cy - dy;
    const adz = az - dz;
    const bdz = bz - dz;
    const cdz = cz - dz;

    return adx * (bdy * cdz - bdz * cdy) +
        bdx * (cdy * adz - cdz * ady) +
        cdx * (ady * bdz - adz * bdy);
}

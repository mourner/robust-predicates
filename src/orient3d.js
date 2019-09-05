import {
    epsilon, splitter, resulterrbound, estimate, vec,
    expansion_sum as sum, scale_expansion as scale
} from './util.js';

const o3derrboundA = (7 + 56 * epsilon) * epsilon;
const o3derrboundB = (3 + 28 * epsilon) * epsilon;
const o3derrboundC = (26 + 288 * epsilon) * epsilon * epsilon;

const bc = vec(4);
const ca = vec(4);
const ab = vec(4);
const adet = vec(8);
const bdet = vec(8);
const cdet = vec(8);
const abdet = vec(16);
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
const v = vec(12);
const w = vec(16);

let fin = vec(192);
let fin2 = vec(192);

function orient3dadapt(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz, permanent) {
    let finlen, tmp;
    let adxtail, bdxtail, cdxtail;
    let adytail, bdytail, cdytail;
    let adztail, bdztail, cdztail;
    let at_blen, at_clen, bt_clen, bt_alen, ct_alen, ct_blen;
    let vlength, wlength;
    let negate;
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
    const alen = scale(4, bc, adz, adet);

    $Cross_Product(cdx, cdy, adx, ady, ca);
    const blen = scale(4, ca, bdz, bdet);

    $Cross_Product(adx, ady, bdx, bdy, ab);
    const clen = scale(4, ab, cdz, cdet);

    const ablen = sum(alen, adet, blen, bdet, abdet);
    finlen = sum(ablen, abdet, clen, cdet, fin);

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

    if (adxtail === 0) {
        if (adytail === 0) {
            at_b[0] = 0;
            at_blen = 1;
            at_c[0] = 0;
            at_clen = 1;
        } else {
            negate = -adytail;
            $Two_Product(negate, bdx, s1, at_b[0]);
            at_b[1] = s1;
            at_blen = 2;
            $Two_Product(adytail, cdx, s1, at_c[0]);
            at_c[1] = s1;
            at_clen = 2;
        }
    } else {
        if (adytail === 0) {
            $Two_Product(adxtail, bdy, s1, at_b[0]);
            at_b[1] = s1;
            at_blen = 2;
            negate = -adxtail;
            $Two_Product(negate, cdy, s1, at_c[0]);
            at_c[1] = s1;
            at_clen = 2;
        } else {
            $Cross_Product(adxtail, bdx, adytail, bdy, at_b);
            at_blen = 4;
            $Cross_Product(adytail, cdy, adxtail, cdx, at_c);
            at_clen = 4;
        }
    }
    if (bdxtail === 0) {
        if (bdytail === 0) {
            bt_c[0] = 0;
            bt_clen = 1;
            bt_a[0] = 0;
            bt_alen = 1;
        } else {
            negate = -bdytail;
            $Two_Product(negate, cdx, s1, bt_c[0]);
            bt_c[1] = s1;
            bt_clen = 2;
            $Two_Product(bdytail, adx, s1, bt_a[0]);
            bt_a[1] = s1;
            bt_alen = 2;
        }
    } else {
        if (bdytail === 0) {
            $Two_Product(bdxtail, cdy, s1, bt_c[0]);
            bt_c[1] = s1;
            bt_clen = 2;
            negate = -bdxtail;
            $Two_Product(negate, ady, s1, bt_a[0]);
            bt_a[1] = s1;
            bt_alen = 2;
        } else {
            $Cross_Product(bdxtail, cdx, bdytail, cdy, bt_c);
            bt_clen = 4;
            $Cross_Product(bdytail, ady, bdxtail, adx, bt_a);
            bt_alen = 4;
        }
    }
    if (cdxtail === 0) {
        if (cdytail === 0) {
            ct_a[0] = 0;
            ct_alen = 1;
            ct_b[0] = 0;
            ct_blen = 1;
        } else {
            negate = -cdytail;
            $Two_Product(negate, adx, s1, ct_a[0]);
            ct_a[1] = s1;
            ct_alen = 2;
            $Two_Product(cdytail, bdx, s1, ct_b[0]);
            ct_b[1] = s1;
            ct_blen = 2;
        }
    } else {
        if (cdytail === 0) {
            $Two_Product(cdxtail, ady, s1, ct_a[0]);
            ct_a[1] = s1;
            ct_alen = 2;
            negate = -cdxtail;
            $Two_Product(negate, bdy, s1, ct_b[0]);
            ct_b[1] = s1;
            ct_blen = 2;
        } else {
            $Cross_Product(cdxtail, adx, cdytail, ady, ct_a);
            ct_alen = 4;
            $Cross_Product(cdytail, bdy, cdxtail, bdx, ct_b);
            ct_blen = 4;
        }
    }

    const bctlen = sum(bt_clen, bt_c, ct_blen, ct_b, bct);
    wlength = scale(bctlen, bct, adz, w);
    finlen = sum(finlen, fin, wlength, w, fin2);
    tmp = fin; fin = fin2; fin2 = tmp;

    const catlen = sum(ct_alen, ct_a, at_clen, at_c, cat);
    wlength = scale(catlen, cat, bdz, w);
    finlen = sum(finlen, fin, wlength, w, fin2);
    tmp = fin; fin = fin2; fin2 = tmp;

    const abtlen = sum(at_blen, at_b, bt_alen, bt_a, abt);
    wlength = scale(abtlen, abt, cdz, w);
    finlen = sum(finlen, fin, wlength, w, fin2);
    tmp = fin; fin = fin2; fin2 = tmp;

    if (adztail !== 0) {
        vlength = scale(4, bc, adztail, v);
        finlen = sum(finlen, fin, vlength, v, fin2);
        tmp = fin; fin = fin2; fin2 = tmp;
    }
    if (bdztail !== 0) {
        vlength = scale(4, ca, bdztail, v);
        finlen = sum(finlen, fin, vlength, v, fin2);
        tmp = fin; fin = fin2; fin2 = tmp;
    }
    if (cdztail !== 0) {
        vlength = scale(4, ab, cdztail, v);
        finlen = sum(finlen, fin, vlength, v, fin2);
        tmp = fin; fin = fin2; fin2 = tmp;
    }

    if (adxtail !== 0) {
        if (bdytail !== 0) {
            $Two_Product(adxtail, bdytail, s1, s0);
            $Two_One_Product(s1, s0, cdz, u);
            finlen = sum(finlen, fin, 4, u, fin2);
            tmp = fin; fin = fin2; fin2 = tmp;
            if (cdztail !== 0) {
                $Two_One_Product(s1, s0, cdztail, u);
                finlen = sum(finlen, fin, 4, u, fin2);
                tmp = fin; fin = fin2; fin2 = tmp;
            }
        }
        if (cdytail !== 0) {
            negate = -adxtail;
            $Two_Product(negate, cdytail, s1, s0);
            $Two_One_Product(s1, s0, bdz, u);
            finlen = sum(finlen, fin, 4, u, fin2);
            tmp = fin; fin = fin2; fin2 = tmp;
            if (bdztail !== 0) {
                $Two_One_Product(s1, s0, bdztail, u);
                finlen = sum(finlen, fin, 4, u, fin2);
                tmp = fin; fin = fin2; fin2 = tmp;
            }
        }
    }
    if (bdxtail !== 0) {
        if (cdytail !== 0) {
            $Two_Product(bdxtail, cdytail, s1, s0);
            $Two_One_Product(s1, s0, adz, u);
            finlen = sum(finlen, fin, 4, u, fin2);
            tmp = fin; fin = fin2; fin2 = tmp;
            if (adztail !== 0) {
                $Two_One_Product(s1, s0, adztail, u);
                finlen = sum(finlen, fin, 4, u, fin2);
                tmp = fin; fin = fin2; fin2 = tmp;
            }
        }
        if (adytail !== 0) {
            negate = -bdxtail;
            $Two_Product(negate, adytail, s1, s0);
            $Two_One_Product(s1, s0, cdz, u);
            finlen = sum(finlen, fin, 4, u, fin2);
            tmp = fin; fin = fin2; fin2 = tmp;
            if (cdztail !== 0) {
                $Two_One_Product(s1, s0, cdztail, u);
                finlen = sum(finlen, fin, 4, u, fin2);
                tmp = fin; fin = fin2; fin2 = tmp;
            }
        }
    }
    if (cdxtail !== 0) {
        if (adytail !== 0) {
            $Two_Product(cdxtail, adytail, s1, s0);
            $Two_One_Product(s1, s0, bdz, u);
            finlen = sum(finlen, fin, 4, u, fin2);
            tmp = fin; fin = fin2; fin2 = tmp;
            if (bdztail !== 0) {
                $Two_One_Product(s1, s0, bdztail, u);
                finlen = sum(finlen, fin, 4, u, fin2);
                tmp = fin; fin = fin2; fin2 = tmp;
            }
        }
        if (bdytail !== 0) {
            negate = -cdxtail;
            $Two_Product(negate, bdytail, s1, s0);
            $Two_One_Product(s1, s0, adz, u);
            finlen = sum(finlen, fin, 4, u, fin2);
            tmp = fin; fin = fin2; fin2 = tmp;
            if (adztail !== 0) {
                $Two_One_Product(s1, s0, adztail, u);
                finlen = sum(finlen, fin, 4, u, fin2);
                tmp = fin; fin = fin2; fin2 = tmp;
            }
        }
    }

    if (adztail !== 0) {
        wlength = scale(bctlen, bct, adztail, w);
        finlen = sum(finlen, fin, wlength, w, fin2);
        tmp = fin; fin = fin2; fin2 = tmp;
    }
    if (bdztail !== 0) {
        wlength = scale(catlen, cat, bdztail, w);
        finlen = sum(finlen, fin, wlength, w, fin2);
        tmp = fin; fin = fin2; fin2 = tmp;
    }
    if (cdztail !== 0) {
        wlength = scale(abtlen, abt, cdztail, w);
        finlen = sum(finlen, fin, wlength, w, fin2);
        tmp = fin; fin = fin2; fin2 = tmp;
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

import {
    epsilon, splitter, resulterrbound, estimate, vec,
    fast_expansion_sum_zeroelim, scale_expansion_zeroelim
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
const fin1 = vec(192);
const fin2 = vec(192);
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

function orient3dadapt(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz, permanent) {
    let finnow, finother, finswap;
    let finlength;
    let adxtail, bdxtail, cdxtail;
    let adytail, bdytail, cdytail;
    let adztail, bdztail, cdztail;
    let at_blarge, at_clarge;
    let bt_clarge, bt_alarge;
    let ct_alarge, ct_blarge;
    let at_blen, at_clen, bt_clen, bt_alen, ct_alen, ct_blen;
    let bdxt_cdyt1, cdxt_bdyt1, cdxt_adyt1;
    let adxt_cdyt1, adxt_bdyt1, bdxt_adyt1;
    let bdxt_cdyt0, cdxt_bdyt0, cdxt_adyt0;
    let adxt_cdyt0, adxt_bdyt0, bdxt_adyt0;
    let s1, s0, t1, t0, u3;
    let vlength, wlength;
    let negate;

    let bvirt, c, ahi, alo, bhi, blo, _i, _j, _k, _0;

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
    const alen = scale_expansion_zeroelim(4, bc, adz, adet);

    $Cross_Product(cdx, cdy, adx, ady, ca);
    const blen = scale_expansion_zeroelim(4, ca, bdz, bdet);

    $Cross_Product(adx, ady, bdx, bdy, ab);
    const clen = scale_expansion_zeroelim(4, ab, cdz, cdet);

    const ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, abdet);
    finlength = fast_expansion_sum_zeroelim(ablen, abdet, clen, cdet, fin1);

    let det = estimate(finlength, fin1);
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

    finnow = fin1;
    finother = fin2;

    if (adxtail === 0) {
        if (adytail === 0) {
            at_b[0] = 0;
            at_blen = 1;
            at_c[0] = 0;
            at_clen = 1;
        } else {
            negate = -adytail;
            $Two_Product(negate, bdx, at_blarge, at_b[0]);
            at_b[1] = at_blarge;
            at_blen = 2;
            $Two_Product(adytail, cdx, at_clarge, at_c[0]);
            at_c[1] = at_clarge;
            at_clen = 2;
        }
    } else {
        if (adytail === 0) {
            $Two_Product(adxtail, bdy, at_blarge, at_b[0]);
            at_b[1] = at_blarge;
            at_blen = 2;
            negate = -adxtail;
            $Two_Product(negate, cdy, at_clarge, at_c[0]);
            at_c[1] = at_clarge;
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
            $Two_Product(negate, cdx, bt_clarge, bt_c[0]);
            bt_c[1] = bt_clarge;
            bt_clen = 2;
            $Two_Product(bdytail, adx, bt_alarge, bt_a[0]);
            bt_a[1] = bt_alarge;
            bt_alen = 2;
        }
    } else {
        if (bdytail === 0) {
            $Two_Product(bdxtail, cdy, bt_clarge, bt_c[0]);
            bt_c[1] = bt_clarge;
            bt_clen = 2;
            negate = -bdxtail;
            $Two_Product(negate, ady, bt_alarge, bt_a[0]);
            bt_a[1] = bt_alarge;
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
            $Two_Product(negate, adx, ct_alarge, ct_a[0]);
            ct_a[1] = ct_alarge;
            ct_alen = 2;
            $Two_Product(cdytail, bdx, ct_blarge, ct_b[0]);
            ct_b[1] = ct_blarge;
            ct_blen = 2;
        }
    } else {
        if (cdytail === 0) {
            $Two_Product(cdxtail, ady, ct_alarge, ct_a[0]);
            ct_a[1] = ct_alarge;
            ct_alen = 2;
            negate = -cdxtail;
            $Two_Product(negate, bdy, ct_blarge, ct_b[0]);
            ct_b[1] = ct_blarge;
            ct_blen = 2;
        } else {
            $Cross_Product(cdxtail, adx, cdytail, ady, ct_a);
            ct_alen = 4;
            $Cross_Product(cdytail, bdy, cdxtail, bdx, ct_b);
            ct_blen = 4;
        }
    }

    const bctlen = fast_expansion_sum_zeroelim(bt_clen, bt_c, ct_blen, ct_b, bct);
    wlength = scale_expansion_zeroelim(bctlen, bct, adz, w);
    finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w, finother);
    finswap = finnow; finnow = finother; finother = finswap;

    const catlen = fast_expansion_sum_zeroelim(ct_alen, ct_a, at_clen, at_c, cat);
    wlength = scale_expansion_zeroelim(catlen, cat, bdz, w);
    finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w, finother);
    finswap = finnow; finnow = finother; finother = finswap;

    const abtlen = fast_expansion_sum_zeroelim(at_blen, at_b, bt_alen, bt_a, abt);
    wlength = scale_expansion_zeroelim(abtlen, abt, cdz, w);
    finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w, finother);
    finswap = finnow; finnow = finother; finother = finswap;

    if (adztail !== 0) {
        vlength = scale_expansion_zeroelim(4, bc, adztail, v);
        finlength = fast_expansion_sum_zeroelim(finlength, finnow, vlength, v, finother);
        finswap = finnow; finnow = finother; finother = finswap;
    }
    if (bdztail !== 0) {
        vlength = scale_expansion_zeroelim(4, ca, bdztail, v);
        finlength = fast_expansion_sum_zeroelim(finlength, finnow, vlength, v, finother);
        finswap = finnow; finnow = finother; finother = finswap;
    }
    if (cdztail !== 0) {
        vlength = scale_expansion_zeroelim(4, ab, cdztail, v);
        finlength = fast_expansion_sum_zeroelim(finlength, finnow, vlength, v, finother);
        finswap = finnow; finnow = finother; finother = finswap;
    }

    if (adxtail !== 0) {
        if (bdytail !== 0) {
            $Two_Product(adxtail, bdytail, adxt_bdyt1, adxt_bdyt0);
            $Two_One_Product(adxt_bdyt1, adxt_bdyt0, cdz, u3, u[2], u[1], u[0]);
            u[3] = u3;
            finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u, finother);
            finswap = finnow; finnow = finother; finother = finswap;
            if (cdztail !== 0) {
                $Two_One_Product(adxt_bdyt1, adxt_bdyt0, cdztail, u3, u[2], u[1], u[0]);
                u[3] = u3;
                finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u, finother);
                finswap = finnow; finnow = finother; finother = finswap;
            }
        }
        if (cdytail !== 0) {
            negate = -adxtail;
            $Two_Product(negate, cdytail, adxt_cdyt1, adxt_cdyt0);
            $Two_One_Product(adxt_cdyt1, adxt_cdyt0, bdz, u3, u[2], u[1], u[0]);
            u[3] = u3;
            finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u, finother);
            finswap = finnow; finnow = finother; finother = finswap;
            if (bdztail !== 0) {
                $Two_One_Product(adxt_cdyt1, adxt_cdyt0, bdztail, u3, u[2], u[1], u[0]);
                u[3] = u3;
                finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u, finother);
                finswap = finnow; finnow = finother; finother = finswap;
            }
        }
    }
    if (bdxtail !== 0) {
        if (cdytail !== 0) {
            $Two_Product(bdxtail, cdytail, bdxt_cdyt1, bdxt_cdyt0);
            $Two_One_Product(bdxt_cdyt1, bdxt_cdyt0, adz, u3, u[2], u[1], u[0]);
            u[3] = u3;
            finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u, finother);
            finswap = finnow; finnow = finother; finother = finswap;
            if (adztail !== 0) {
                $Two_One_Product(bdxt_cdyt1, bdxt_cdyt0, adztail, u3, u[2], u[1], u[0]);
                u[3] = u3;
                finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u, finother);
                finswap = finnow; finnow = finother; finother = finswap;
            }
        }
        if (adytail !== 0) {
            negate = -bdxtail;
            $Two_Product(negate, adytail, bdxt_adyt1, bdxt_adyt0);
            $Two_One_Product(bdxt_adyt1, bdxt_adyt0, cdz, u3, u[2], u[1], u[0]);
            u[3] = u3;
            finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u, finother);
            finswap = finnow; finnow = finother; finother = finswap;
            if (cdztail !== 0) {
                $Two_One_Product(bdxt_adyt1, bdxt_adyt0, cdztail, u3, u[2], u[1], u[0]);
                u[3] = u3;
                finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u, finother);
                finswap = finnow; finnow = finother; finother = finswap;
            }
        }
    }
    if (cdxtail !== 0) {
        if (adytail !== 0) {
            $Two_Product(cdxtail, adytail, cdxt_adyt1, cdxt_adyt0);
            $Two_One_Product(cdxt_adyt1, cdxt_adyt0, bdz, u3, u[2], u[1], u[0]);
            u[3] = u3;
            finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u, finother);
            finswap = finnow; finnow = finother; finother = finswap;
            if (bdztail !== 0) {
                $Two_One_Product(cdxt_adyt1, cdxt_adyt0, bdztail, u3, u[2], u[1], u[0]);
                u[3] = u3;
                finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u, finother);
                finswap = finnow; finnow = finother; finother = finswap;
            }
        }
        if (bdytail !== 0) {
            negate = -cdxtail;
            $Two_Product(negate, bdytail, cdxt_bdyt1, cdxt_bdyt0);
            $Two_One_Product(cdxt_bdyt1, cdxt_bdyt0, adz, u3, u[2], u[1], u[0]);
            u[3] = u3;
            finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u, finother);
            finswap = finnow; finnow = finother; finother = finswap;
            if (adztail !== 0) {
                $Two_One_Product(cdxt_bdyt1, cdxt_bdyt0, adztail, u3, u[2], u[1], u[0]);
                u[3] = u3;
                finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u, finother);
                finswap = finnow; finnow = finother; finother = finswap;
            }
        }
    }

    if (adztail !== 0) {
        wlength = scale_expansion_zeroelim(bctlen, bct, adztail, w);
        finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w, finother);
        finswap = finnow; finnow = finother; finother = finswap;
    }
    if (bdztail !== 0) {
        wlength = scale_expansion_zeroelim(catlen, cat, bdztail, w);
        finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w, finother);
        finswap = finnow; finnow = finother; finother = finswap;
    }
    if (cdztail !== 0) {
        wlength = scale_expansion_zeroelim(abtlen, abt, cdztail, w);
        finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w, finother);
        finswap = finnow; finnow = finother; finother = finswap;
    }

    return finnow[finlength - 1];
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


import {
    epsilon, splitter, resulterrbound, estimate,
    fast_expansion_sum_zeroelim, scale_expansion_zeroelim
} from './util.js';

const o3derrboundA = (7.0 + 56.0 * epsilon) * epsilon;
const o3derrboundB = (3.0 + 28.0 * epsilon) * epsilon;
const o3derrboundC = (26.0 + 288.0 * epsilon) * epsilon * epsilon;

const bc = new Float64Array(4);
const ca = new Float64Array(4);
const ab = new Float64Array(4);
const adet = new Float64Array(8);
const bdet = new Float64Array(8);
const cdet = new Float64Array(8);
const abdet = new Float64Array(16);
const fin1 = new Float64Array(192);
const fin2 = new Float64Array(192);
const at_b = new Float64Array(4);
const at_c = new Float64Array(4);
const bt_c = new Float64Array(4);
const bt_a = new Float64Array(4);
const ct_a = new Float64Array(4);
const ct_b = new Float64Array(4);
const bct = new Float64Array(8);
const cat = new Float64Array(8);
const abt = new Float64Array(8);
const u = new Float64Array(4);
const v = new Float64Array(12);
const w = new Float64Array(16);

function orient3dadapt(pax, pay, paz, pbx, pby, pbz, pcx, pcy, pcz, pdx, pdy, pdz, permanent) {
    let bdxcdy1, cdxbdy1, cdxady1, adxcdy1, adxbdy1, bdxady1;
    let bdxcdy0, cdxbdy0, cdxady0, adxcdy0, adxbdy0, bdxady0;
    let bc3, ca3, ab3;
    let alen, blen, clen;
    let ablen;
    let finnow, finother, finswap;
    let finlength;

    let adxtail, bdxtail, cdxtail;
    let adytail, bdytail, cdytail;
    let adztail, bdztail, cdztail;
    let at_blarge, at_clarge;
    let bt_clarge, bt_alarge;
    let ct_alarge, ct_blarge;
    let at_blen, at_clen, bt_clen, bt_alen, ct_alen, ct_blen;
    let bdxt_cdy1, cdxt_bdy1, cdxt_ady1;
    let adxt_cdy1, adxt_bdy1, bdxt_ady1;
    let bdxt_cdy0, cdxt_bdy0, cdxt_ady0;
    let adxt_cdy0, adxt_bdy0, bdxt_ady0;
    let bdyt_cdx1, cdyt_bdx1, cdyt_adx1;
    let adyt_cdx1, adyt_bdx1, bdyt_adx1;
    let bdyt_cdx0, cdyt_bdx0, cdyt_adx0;
    let adyt_cdx0, adyt_bdx0, bdyt_adx0;
    let bctlen, catlen, abtlen;
    let bdxt_cdyt1, cdxt_bdyt1, cdxt_adyt1;
    let adxt_cdyt1, adxt_bdyt1, bdxt_adyt1;
    let bdxt_cdyt0, cdxt_bdyt0, cdxt_adyt0;
    let adxt_cdyt0, adxt_bdyt0, bdxt_adyt0;
    let u3;
    let vlength, wlength;
    let negate;

    let bvirt, c, ahi, alo, bhi, blo, _i, _j, _k, _0;

    const adx = pax - pdx;
    const bdx = pbx - pdx;
    const cdx = pcx - pdx;
    const ady = pay - pdy;
    const bdy = pby - pdy;
    const cdy = pcy - pdy;
    const adz = paz - pdz;
    const bdz = pbz - pdz;
    const cdz = pcz - pdz;

    $Two_Product(bdx, cdy, bdxcdy1, bdxcdy0);
    $Two_Product(cdx, bdy, cdxbdy1, cdxbdy0);
    $Two_Two_Diff(bdxcdy1, bdxcdy0, cdxbdy1, cdxbdy0, bc3, bc[2], bc[1], bc[0]);
    bc[3] = bc3;
    alen = scale_expansion_zeroelim(4, bc, adz, adet);

    $Two_Product(cdx, ady, cdxady1, cdxady0);
    $Two_Product(adx, cdy, adxcdy1, adxcdy0);
    $Two_Two_Diff(cdxady1, cdxady0, adxcdy1, adxcdy0, ca3, ca[2], ca[1], ca[0]);
    ca[3] = ca3;
    blen = scale_expansion_zeroelim(4, ca, bdz, bdet);

    $Two_Product(adx, bdy, adxbdy1, adxbdy0);
    $Two_Product(bdx, ady, bdxady1, bdxady0);
    $Two_Two_Diff(adxbdy1, adxbdy0, bdxady1, bdxady0, ab3, ab[2], ab[1], ab[0]);
    ab[3] = ab3;
    clen = scale_expansion_zeroelim(4, ab, cdz, cdet);

    ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, abdet);
    finlength = fast_expansion_sum_zeroelim(ablen, abdet, clen, cdet, fin1);

    let det = estimate(finlength, fin1);
    let errbound = o3derrboundB * permanent;
    if ((det >= errbound) || (-det >= errbound)) {
        return det;
    }

    $Two_Diff_Tail(pax, pdx, adx, adxtail);
    $Two_Diff_Tail(pbx, pdx, bdx, bdxtail);
    $Two_Diff_Tail(pcx, pdx, cdx, cdxtail);
    $Two_Diff_Tail(pay, pdy, ady, adytail);
    $Two_Diff_Tail(pby, pdy, bdy, bdytail);
    $Two_Diff_Tail(pcy, pdy, cdy, cdytail);
    $Two_Diff_Tail(paz, pdz, adz, adztail);
    $Two_Diff_Tail(pbz, pdz, bdz, bdztail);
    $Two_Diff_Tail(pcz, pdz, cdz, cdztail);

    if (
        (adxtail === 0.0) && (bdxtail === 0.0) && (cdxtail === 0.0) &&
        (adytail === 0.0) && (bdytail === 0.0) && (cdytail === 0.0) &&
        (adztail === 0.0) && (bdztail === 0.0) && (cdztail === 0.0)) {
        return det;
    }

    errbound = o3derrboundC * permanent + resulterrbound * Math.abs(det);
    det +=
        adz * (bdx * cdytail + cdy * bdxtail - (bdy * cdxtail + cdx * bdytail)) + adztail * (bdx * cdy - bdy * cdx) +
        bdz * (cdx * adytail + ady * cdxtail - (cdy * adxtail + adx * cdytail)) + bdztail * (cdx * ady - cdy * adx) +
        cdz * (adx * bdytail + bdy * adxtail - (ady * bdxtail + bdx * adytail)) + cdztail * (adx * bdy - ady * bdx);
    if ((det >= errbound) || (-det >= errbound)) {
        return det;
    }

    finnow = fin1;
    finother = fin2;

    if (adxtail === 0.0) {
        if (adytail === 0.0) {
            at_b[0] = 0.0;
            at_blen = 1;
            at_c[0] = 0.0;
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
        if (adytail === 0.0) {
            $Two_Product(adxtail, bdy, at_blarge, at_b[0]);
            at_b[1] = at_blarge;
            at_blen = 2;
            negate = -adxtail;
            $Two_Product(negate, cdy, at_clarge, at_c[0]);
            at_c[1] = at_clarge;
            at_clen = 2;
        } else {
            $Two_Product(adxtail, bdy, adxt_bdy1, adxt_bdy0);
            $Two_Product(adytail, bdx, adyt_bdx1, adyt_bdx0);
            $Two_Two_Diff(adxt_bdy1, adxt_bdy0, adyt_bdx1, adyt_bdx0, at_blarge, at_b[2], at_b[1], at_b[0]);
            at_b[3] = at_blarge;
            at_blen = 4;
            $Two_Product(adytail, cdx, adyt_cdx1, adyt_cdx0);
            $Two_Product(adxtail, cdy, adxt_cdy1, adxt_cdy0);
            $Two_Two_Diff(adyt_cdx1, adyt_cdx0, adxt_cdy1, adxt_cdy0, at_clarge, at_c[2], at_c[1], at_c[0]);
            at_c[3] = at_clarge;
            at_clen = 4;
        }
    }
    if (bdxtail === 0.0) {
        if (bdytail === 0.0) {
            bt_c[0] = 0.0;
            bt_clen = 1;
            bt_a[0] = 0.0;
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
        if (bdytail === 0.0) {
            $Two_Product(bdxtail, cdy, bt_clarge, bt_c[0]);
            bt_c[1] = bt_clarge;
            bt_clen = 2;
            negate = -bdxtail;
            $Two_Product(negate, ady, bt_alarge, bt_a[0]);
            bt_a[1] = bt_alarge;
            bt_alen = 2;
        } else {
            $Two_Product(bdxtail, cdy, bdxt_cdy1, bdxt_cdy0);
            $Two_Product(bdytail, cdx, bdyt_cdx1, bdyt_cdx0);
            $Two_Two_Diff(bdxt_cdy1, bdxt_cdy0, bdyt_cdx1, bdyt_cdx0, bt_clarge, bt_c[2], bt_c[1], bt_c[0]);
            bt_c[3] = bt_clarge;
            bt_clen = 4;
            $Two_Product(bdytail, adx, bdyt_adx1, bdyt_adx0);
            $Two_Product(bdxtail, ady, bdxt_ady1, bdxt_ady0);
            $Two_Two_Diff(bdyt_adx1, bdyt_adx0, bdxt_ady1, bdxt_ady0, bt_alarge, bt_a[2], bt_a[1], bt_a[0]);
            bt_a[3] = bt_alarge;
            bt_alen = 4;
        }
    }
    if (cdxtail === 0.0) {
        if (cdytail === 0.0) {
            ct_a[0] = 0.0;
            ct_alen = 1;
            ct_b[0] = 0.0;
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
        if (cdytail === 0.0) {
            $Two_Product(cdxtail, ady, ct_alarge, ct_a[0]);
            ct_a[1] = ct_alarge;
            ct_alen = 2;
            negate = -cdxtail;
            $Two_Product(negate, bdy, ct_blarge, ct_b[0]);
            ct_b[1] = ct_blarge;
            ct_blen = 2;
        } else {
            $Two_Product(cdxtail, ady, cdxt_ady1, cdxt_ady0);
            $Two_Product(cdytail, adx, cdyt_adx1, cdyt_adx0);
            $Two_Two_Diff(cdxt_ady1, cdxt_ady0, cdyt_adx1, cdyt_adx0, ct_alarge, ct_a[2], ct_a[1], ct_a[0]);
            ct_a[3] = ct_alarge;
            ct_alen = 4;
            $Two_Product(cdytail, bdx, cdyt_bdx1, cdyt_bdx0);
            $Two_Product(cdxtail, bdy, cdxt_bdy1, cdxt_bdy0);
            $Two_Two_Diff(cdyt_bdx1, cdyt_bdx0, cdxt_bdy1, cdxt_bdy0, ct_blarge, ct_b[2], ct_b[1], ct_b[0]);
            ct_b[3] = ct_blarge;
            ct_blen = 4;
        }
    }

    bctlen = fast_expansion_sum_zeroelim(bt_clen, bt_c, ct_blen, ct_b, bct);
    wlength = scale_expansion_zeroelim(bctlen, bct, adz, w);
    finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w, finother);
    finswap = finnow; finnow = finother; finother = finswap;

    catlen = fast_expansion_sum_zeroelim(ct_alen, ct_a, at_clen, at_c, cat);
    wlength = scale_expansion_zeroelim(catlen, cat, bdz, w);
    finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w, finother);
    finswap = finnow; finnow = finother; finother = finswap;

    abtlen = fast_expansion_sum_zeroelim(at_blen, at_b, bt_alen, bt_a, abt);
    wlength = scale_expansion_zeroelim(abtlen, abt, cdz, w);
    finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w, finother);
    finswap = finnow; finnow = finother; finother = finswap;

    if (adztail !== 0.0) {
        vlength = scale_expansion_zeroelim(4, bc, adztail, v);
        finlength = fast_expansion_sum_zeroelim(finlength, finnow, vlength, v, finother);
        finswap = finnow; finnow = finother; finother = finswap;
    }
    if (bdztail !== 0.0) {
        vlength = scale_expansion_zeroelim(4, ca, bdztail, v);
        finlength = fast_expansion_sum_zeroelim(finlength, finnow, vlength, v, finother);
        finswap = finnow; finnow = finother; finother = finswap;
    }
    if (cdztail !== 0.0) {
        vlength = scale_expansion_zeroelim(4, ab, cdztail, v);
        finlength = fast_expansion_sum_zeroelim(finlength, finnow, vlength, v, finother);
        finswap = finnow; finnow = finother; finother = finswap;
    }

    if (adxtail !== 0.0) {
        if (bdytail !== 0.0) {
            $Two_Product(adxtail, bdytail, adxt_bdyt1, adxt_bdyt0);
            $Two_One_Product(adxt_bdyt1, adxt_bdyt0, cdz, u3, u[2], u[1], u[0]);
            u[3] = u3;
            finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u, finother);
            finswap = finnow; finnow = finother; finother = finswap;
            if (cdztail !== 0.0) {
                $Two_One_Product(adxt_bdyt1, adxt_bdyt0, cdztail, u3, u[2], u[1], u[0]);
                u[3] = u3;
                finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u, finother);
                finswap = finnow; finnow = finother; finother = finswap;
            }
        }
        if (cdytail !== 0.0) {
            negate = -adxtail;
            $Two_Product(negate, cdytail, adxt_cdyt1, adxt_cdyt0);
            $Two_One_Product(adxt_cdyt1, adxt_cdyt0, bdz, u3, u[2], u[1], u[0]);
            u[3] = u3;
            finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u, finother);
            finswap = finnow; finnow = finother; finother = finswap;
            if (bdztail !== 0.0) {
                $Two_One_Product(adxt_cdyt1, adxt_cdyt0, bdztail, u3, u[2], u[1], u[0]);
                u[3] = u3;
                finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u, finother);
                finswap = finnow; finnow = finother; finother = finswap;
            }
        }
    }
    if (bdxtail !== 0.0) {
        if (cdytail !== 0.0) {
            $Two_Product(bdxtail, cdytail, bdxt_cdyt1, bdxt_cdyt0);
            $Two_One_Product(bdxt_cdyt1, bdxt_cdyt0, adz, u3, u[2], u[1], u[0]);
            u[3] = u3;
            finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u, finother);
            finswap = finnow; finnow = finother; finother = finswap;
            if (adztail !== 0.0) {
                $Two_One_Product(bdxt_cdyt1, bdxt_cdyt0, adztail, u3, u[2], u[1], u[0]);
                u[3] = u3;
                finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u, finother);
                finswap = finnow; finnow = finother; finother = finswap;
            }
        }
        if (adytail !== 0.0) {
            negate = -bdxtail;
            $Two_Product(negate, adytail, bdxt_adyt1, bdxt_adyt0);
            $Two_One_Product(bdxt_adyt1, bdxt_adyt0, cdz, u3, u[2], u[1], u[0]);
            u[3] = u3;
            finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u, finother);
            finswap = finnow; finnow = finother; finother = finswap;
            if (cdztail !== 0.0) {
                $Two_One_Product(bdxt_adyt1, bdxt_adyt0, cdztail, u3, u[2], u[1], u[0]);
                u[3] = u3;
                finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u, finother);
                finswap = finnow; finnow = finother; finother = finswap;
            }
        }
    }
    if (cdxtail !== 0.0) {
        if (adytail !== 0.0) {
            $Two_Product(cdxtail, adytail, cdxt_adyt1, cdxt_adyt0);
            $Two_One_Product(cdxt_adyt1, cdxt_adyt0, bdz, u3, u[2], u[1], u[0]);
            u[3] = u3;
            finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u, finother);
            finswap = finnow; finnow = finother; finother = finswap;
            if (bdztail !== 0.0) {
                $Two_One_Product(cdxt_adyt1, cdxt_adyt0, bdztail, u3, u[2], u[1], u[0]);
                u[3] = u3;
                finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u, finother);
                finswap = finnow; finnow = finother; finother = finswap;
            }
        }
        if (bdytail !== 0.0) {
            negate = -cdxtail;
            $Two_Product(negate, bdytail, cdxt_bdyt1, cdxt_bdyt0);
            $Two_One_Product(cdxt_bdyt1, cdxt_bdyt0, adz, u3, u[2], u[1], u[0]);
            u[3] = u3;
            finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u, finother);
            finswap = finnow; finnow = finother; finother = finswap;
            if (adztail !== 0.0) {
                $Two_One_Product(cdxt_bdyt1, cdxt_bdyt0, adztail, u3, u[2], u[1], u[0]);
                u[3] = u3;
                finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u, finother);
                finswap = finnow; finnow = finother; finother = finswap;
            }
        }
    }

    if (adztail !== 0.0) {
        wlength = scale_expansion_zeroelim(bctlen, bct, adztail, w);
        finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w, finother);
        finswap = finnow; finnow = finother; finother = finswap;
    }
    if (bdztail !== 0.0) {
        wlength = scale_expansion_zeroelim(catlen, cat, bdztail, w);
        finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w, finother);
        finswap = finnow; finnow = finother; finother = finswap;
    }
    if (cdztail !== 0.0) {
        wlength = scale_expansion_zeroelim(abtlen, abt, cdztail, w);
        finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w, finother);
        finswap = finnow; finnow = finother; finother = finswap;
    }

    return finnow[finlength - 1];
}

export function orient3d(pax, pay, paz, pbx, pby, pbz, pcx, pcy, pcz, pdx, pdy, pdz) {
    const adx = pax - pdx;
    const bdx = pbx - pdx;
    const cdx = pcx - pdx;
    const ady = pay - pdy;
    const bdy = pby - pdy;
    const cdy = pcy - pdy;
    const adz = paz - pdz;
    const bdz = pbz - pdz;
    const cdz = pcz - pdz;

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
    if ((det > errbound) || (-det > errbound)) {
        return det;
    }

    return orient3dadapt(pax, pay, paz, pbx, pby, pbz, pcx, pcy, pcz, pdx, pdy, pdz, permanent);
}

import {
    epsilon, splitter, resulterrbound, estimate, vec,
    expansion_sum, scale_expansion, scale_expansion_twice
} from './util.js';

const iccerrboundA = (10 + 96 * epsilon) * epsilon;
const iccerrboundB = (4 + 48 * epsilon) * epsilon;
const iccerrboundC = (44 + 576 * epsilon) * epsilon * epsilon;

const bc = vec(4);
const ca = vec(4);
const ab = vec(4);
const fin1 = vec(1152);
const fin2 = vec(1152);
const aa = vec(4);
const bb = vec(4);
const cc = vec(4);
const u = vec(4);
const v = vec(4);
const temp8 = vec(8);
const temp16a = vec(16);
const temp16b = vec(16);
const temp16c = vec(16);
const temp32a = vec(32);
const temp32b = vec(32);
const temp48 = vec(48);
const temp64 = vec(64);
const axtbc = vec(8);
const aytbc = vec(8);
const bxtca = vec(8);
const bytca = vec(8);
const cxtab = vec(8);
const cytab = vec(8);
const abt = vec(8);
const bct = vec(8);
const cat = vec(8);
const abtt = vec(4);
const bctt = vec(4);
const catt = vec(4);

function incircleadapt(ax, ay, bx, by, cx, cy, dx, dy, permanent) {
    let finnow, finother, finswap, finlength;
    let adxtail, bdxtail, cdxtail, adytail, bdytail, cdytail;
    let temp8len, temp16alen, temp16blen, temp16clen;
    let temp32alen, temp32blen, temp48len, temp64len;
    let axtbclen, aytbclen, bxtcalen, bytcalen, cxtablen, cytablen;
    let abtlen, bctlen, catlen;
    let abttlen, bcttlen, cattlen;
    let n1, n0;

    let bvirt, c, ahi, alo, bhi, blo, _i, _j, _0, s1, s0, t1, t0, u3;

    const adx = ax - dx;
    const bdx = bx - dx;
    const cdx = cx - dx;
    const ady = ay - dy;
    const bdy = by - dy;
    const cdy = cy - dy;

    $Cross_Product(bdx, bdy, cdx, cdy, bc);
    temp16alen = scale_expansion_twice(4, bc, adx, temp8, adx, temp16a);
    temp16blen = scale_expansion_twice(4, bc, ady, temp8, ady, temp16b);
    temp32alen = expansion_sum(temp16alen, temp16a, temp16blen, temp16b, temp32a);

    $Cross_Product(cdx, cdy, adx, ady, ca);
    temp16alen = scale_expansion_twice(4, ca, bdx, temp8, bdx, temp16a);
    temp16blen = scale_expansion_twice(4, ca, bdy, temp8, bdy, temp16b);
    temp32blen = expansion_sum(temp16alen, temp16a, temp16blen, temp16b, temp32b);

    temp64len = expansion_sum(temp32alen, temp32a, temp32blen, temp32b, temp64);

    $Cross_Product(adx, ady, bdx, bdy, ab);
    temp16alen = scale_expansion_twice(4, ab, cdx, temp8, cdx, temp16a);
    temp16blen = scale_expansion_twice(4, ab, cdy, temp8, cdy, temp16b);
    const clen = expansion_sum(temp16alen, temp16a, temp16blen, temp16b, temp32a);

    finlength = expansion_sum(temp64len, temp64, clen, temp32a, fin1);

    let det = estimate(finlength, fin1);
    let errbound = iccerrboundB * permanent;
    if (det >= errbound || -det >= errbound) {
        return det;
    }

    $Two_Diff_Tail(ax, dx, adx, adxtail);
    $Two_Diff_Tail(ay, dy, ady, adytail);
    $Two_Diff_Tail(bx, dx, bdx, bdxtail);
    $Two_Diff_Tail(by, dy, bdy, bdytail);
    $Two_Diff_Tail(cx, dx, cdx, cdxtail);
    $Two_Diff_Tail(cy, dy, cdy, cdytail);
    if (adxtail === 0 && bdxtail === 0 && cdxtail === 0 && adytail === 0 && bdytail === 0 && cdytail === 0) {
        return det;
    }

    errbound = iccerrboundC * permanent + resulterrbound * Math.abs(det);
    det += ((adx * adx + ady * ady) * ((bdx * cdytail + cdy * bdxtail) - (bdy * cdxtail + cdx * bdytail)) +
        2 * (adx * adxtail + ady * adytail) * (bdx * cdy - bdy * cdx)) +
        ((bdx * bdx + bdy * bdy) * ((cdx * adytail + ady * cdxtail) - (cdy * adxtail + adx * cdytail)) +
        2 * (bdx * bdxtail + bdy * bdytail) * (cdx * ady - cdy * adx)) +
        ((cdx * cdx + cdy * cdy) * ((adx * bdytail + bdy * adxtail) - (ady * bdxtail + bdx * adytail)) +
        2 * (cdx * cdxtail + cdy * cdytail) * (adx * bdy - ady * bdx));

    if (det >= errbound || -det >= errbound) {
        return det;
    }

    finnow = fin1;
    finother = fin2;

    if (bdxtail !== 0 || bdytail !== 0 || cdxtail !== 0 || cdytail !== 0) {
        $Square_Sum(adx, ady, aa);
    }
    if (cdxtail !== 0 || cdytail !== 0 || adxtail !== 0 || adytail !== 0) {
        $Square_Sum(bdx, bdy, bb);
    }
    if (adxtail !== 0 || adytail !== 0 || bdxtail !== 0 || bdytail !== 0) {
        $Square_Sum(cdx, cdy, cc);
    }

    if (adxtail !== 0) {
        axtbclen = scale_expansion(4, bc, adxtail, axtbc);
        temp16alen = scale_expansion(axtbclen, axtbc, 2 * adx, temp16a);
        temp16blen = scale_expansion_twice(4, cc, adxtail, temp8, bdy, temp16b);
        temp16clen = scale_expansion_twice(4, bb, adxtail, temp8, -cdy, temp16c);
        temp32alen = expansion_sum(temp16alen, temp16a, temp16blen, temp16b, temp32a);
        temp48len = expansion_sum(temp16clen, temp16c, temp32alen, temp32a, temp48);
        finlength = expansion_sum(finlength, finnow, temp48len, temp48, finother);
        finswap = finnow; finnow = finother; finother = finswap;
    }
    if (adytail !== 0) {
        aytbclen = scale_expansion(4, bc, adytail, aytbc);
        temp16alen = scale_expansion(aytbclen, aytbc, 2 * ady, temp16a);
        temp16blen = scale_expansion_twice(4, bb, adytail, temp8, cdx, temp16b);
        temp16clen = scale_expansion_twice(4, cc, adytail, temp8, -bdx, temp16c);
        temp32alen = expansion_sum(temp16alen, temp16a, temp16blen, temp16b, temp32a);
        temp48len = expansion_sum(temp16clen, temp16c, temp32alen, temp32a, temp48);
        finlength = expansion_sum(finlength, finnow, temp48len, temp48, finother);
        finswap = finnow; finnow = finother; finother = finswap;
    }
    if (bdxtail !== 0) {
        bxtcalen = scale_expansion(4, ca, bdxtail, bxtca);
        temp16alen = scale_expansion(bxtcalen, bxtca, 2 * bdx, temp16a);
        temp16blen = scale_expansion_twice(4, aa, bdxtail, temp8, cdy, temp16b);
        temp16clen = scale_expansion_twice(4, cc, bdxtail, temp8, -ady, temp16c);
        temp32alen = expansion_sum(temp16alen, temp16a, temp16blen, temp16b, temp32a);
        temp48len = expansion_sum(temp16clen, temp16c, temp32alen, temp32a, temp48);
        finlength = expansion_sum(finlength, finnow, temp48len, temp48, finother);
        finswap = finnow; finnow = finother; finother = finswap;
    }
    if (bdytail !== 0) {
        bytcalen = scale_expansion(4, ca, bdytail, bytca);
        temp16alen = scale_expansion(bytcalen, bytca, 2 * bdy, temp16a);
        temp16blen = scale_expansion_twice(4, cc, bdytail, temp8, adx, temp16b);
        temp16clen = scale_expansion_twice(4, aa, bdytail, temp8, -cdx, temp16c);
        temp32alen = expansion_sum(temp16alen, temp16a, temp16blen, temp16b, temp32a);
        temp48len = expansion_sum(temp16clen, temp16c, temp32alen, temp32a, temp48);
        finlength = expansion_sum(finlength, finnow, temp48len, temp48, finother);
        finswap = finnow; finnow = finother; finother = finswap;
    }
    if (cdxtail !== 0) {
        cxtablen = scale_expansion(4, ab, cdxtail, cxtab);
        temp16alen = scale_expansion(cxtablen, cxtab, 2 * cdx, temp16a);
        temp16blen = scale_expansion_twice(4, bb, cdxtail, temp8, ady, temp16b);
        temp16clen = scale_expansion_twice(4, aa, cdxtail, temp8, -bdy, temp16c);
        temp32alen = expansion_sum(temp16alen, temp16a, temp16blen, temp16b, temp32a);
        temp48len = expansion_sum(temp16clen, temp16c, temp32alen, temp32a, temp48);
        finlength = expansion_sum(finlength, finnow, temp48len, temp48, finother);
        finswap = finnow; finnow = finother; finother = finswap;
    }
    if (cdytail !== 0) {
        cytablen = scale_expansion(4, ab, cdytail, cytab);
        temp16alen = scale_expansion(cytablen, cytab, 2 * cdy, temp16a);
        temp16blen = scale_expansion_twice(4, aa, cdytail, temp8, bdx, temp16b);
        temp16clen = scale_expansion_twice(4, bb, cdytail, temp8, -adx, temp16c);
        temp32alen = expansion_sum(temp16alen, temp16a, temp16blen, temp16b, temp32a);
        temp48len = expansion_sum(temp16clen, temp16c, temp32alen, temp32a, temp48);
        finlength = expansion_sum(finlength, finnow, temp48len, temp48, finother);
        finswap = finnow; finnow = finother; finother = finswap;
    }

    if (adxtail !== 0 || adytail !== 0) {
        if (bdxtail !== 0 || bdytail !== 0 || cdxtail !== 0 || cdytail !== 0) {
            $Two_Product_Sum(bdxtail, cdy, bdx, cdytail, u);
            n1 = -bdy;
            n0 = -bdytail;
            $Two_Product_Sum(cdxtail, n1, cdx, n0, v);
            bctlen = expansion_sum(4, u, 4, v, bct);

            $Cross_Product(bdxtail, bdytail, cdxtail, cdytail, bctt);
            bcttlen = 4;
        } else {
            bct[0] = 0;
            bctlen = 1;
            bctt[0] = 0;
            bcttlen = 1;
        }

        if (adxtail !== 0) {
            temp16alen = scale_expansion(axtbclen, axtbc, adxtail, temp16a);
            temp16blen = scale_expansion(bctlen, bct, adxtail, temp16b);
            temp32alen = scale_expansion(temp16blen, temp16b, 2 * adx, temp32a);
            temp48len = expansion_sum(temp16alen, temp16a, temp32alen, temp32a, temp48);
            finlength = expansion_sum(finlength, finnow, temp48len, temp48, finother);
            finswap = finnow; finnow = finother; finother = finswap;

            if (bdytail !== 0) {
                temp16alen = scale_expansion_twice(4, cc, adxtail, temp8, bdytail, temp16a);
                finlength = expansion_sum(finlength, finnow, temp16alen, temp16a, finother);
                finswap = finnow; finnow = finother; finother = finswap;
            }
            if (cdytail !== 0) {
                temp16alen = scale_expansion_twice(4, bb, -adxtail, temp8, cdytail, temp16a);
                finlength = expansion_sum(finlength, finnow, temp16alen, temp16a, finother);
                finswap = finnow; finnow = finother; finother = finswap;
            }

            temp32alen = scale_expansion(temp16blen, temp16b, adxtail, temp32a);
            temp8len = scale_expansion(bcttlen, bctt, adxtail, temp8);
            temp16alen = scale_expansion(temp8len, temp8, 2 * adx, temp16a);
            temp16blen = scale_expansion(temp8len, temp8, adxtail, temp16b);
            temp32blen = expansion_sum(temp16alen, temp16a, temp16blen, temp16b, temp32b);
            temp64len = expansion_sum(temp32alen, temp32a, temp32blen, temp32b, temp64);
            finlength = expansion_sum(finlength, finnow, temp64len, temp64, finother);
            finswap = finnow; finnow = finother; finother = finswap;
        }

        if (adytail !== 0) {
            temp16alen = scale_expansion(aytbclen, aytbc, adytail, temp16a);
            temp16blen = scale_expansion(bctlen, bct, adytail, temp16b);
            temp32alen = scale_expansion(temp16blen, temp16b, 2 * ady, temp32a);
            temp48len = expansion_sum(temp16alen, temp16a, temp32alen, temp32a, temp48);
            finlength = expansion_sum(finlength, finnow, temp48len, temp48, finother);
            finswap = finnow; finnow = finother; finother = finswap;

            temp32alen = scale_expansion(temp16blen, temp16b, adytail, temp32a);
            temp8len = scale_expansion(bcttlen, bctt, adytail, temp8);
            temp16alen = scale_expansion(temp8len, temp8, 2 * ady, temp16a);
            temp16blen = scale_expansion(temp8len, temp8, adytail, temp16b);
            temp32blen = expansion_sum(temp16alen, temp16a, temp16blen, temp16b, temp32b);
            temp64len = expansion_sum(temp32alen, temp32a, temp32blen, temp32b, temp64);
            finlength = expansion_sum(finlength, finnow, temp64len, temp64, finother);
            finswap = finnow; finnow = finother; finother = finswap;
        }
    }
    if (bdxtail !== 0 || bdytail !== 0) {
        if (cdxtail !== 0 || cdytail !== 0 || adxtail !== 0 || adytail !== 0) {
            $Two_Product_Sum(cdxtail, ady, cdx, adytail, u);
            n1 = -cdy;
            n0 = -cdytail;
            $Two_Product_Sum(adxtail, n1, adx, n0, v);
            catlen = expansion_sum(4, u, 4, v, cat);

            $Cross_Product(cdxtail, cdytail, adxtail, adytail, catt);
            cattlen = 4;
        } else {
            cat[0] = 0;
            catlen = 1;
            catt[0] = 0;
            cattlen = 1;
        }

        if (bdxtail !== 0) {
            temp16alen = scale_expansion(bxtcalen, bxtca, bdxtail, temp16a);
            temp16blen = scale_expansion(catlen, cat, bdxtail, temp16b);
            temp32alen = scale_expansion(temp16blen, temp16b, 2 * bdx, temp32a);
            temp48len = expansion_sum(temp16alen, temp16a, temp32alen, temp32a, temp48);
            finlength = expansion_sum(finlength, finnow, temp48len, temp48, finother);
            finswap = finnow; finnow = finother; finother = finswap;
            if (cdytail !== 0) {
                temp16alen = scale_expansion_twice(4, aa, bdxtail, temp8, cdytail, temp16a);
                finlength = expansion_sum(finlength, finnow, temp16alen, temp16a, finother);
                finswap = finnow; finnow = finother; finother = finswap;
            }
            if (adytail !== 0) {
                temp16alen = scale_expansion_twice(4, cc, -bdxtail, temp8, adytail, temp16a);
                finlength = expansion_sum(finlength, finnow, temp16alen, temp16a, finother);
                finswap = finnow; finnow = finother; finother = finswap;
            }

            temp32alen = scale_expansion(temp16blen, temp16b, bdxtail, temp32a);
            temp8len = scale_expansion(cattlen, catt, bdxtail, temp8);
            temp16alen = scale_expansion(temp8len, temp8, 2 * bdx, temp16a);
            temp16blen = scale_expansion(temp8len, temp8, bdxtail, temp16b);
            temp32blen = expansion_sum(temp16alen, temp16a, temp16blen, temp16b, temp32b);
            temp64len = expansion_sum(temp32alen, temp32a, temp32blen, temp32b, temp64);
            finlength = expansion_sum(finlength, finnow, temp64len, temp64, finother);
            finswap = finnow; finnow = finother; finother = finswap;
        }
        if (bdytail !== 0) {
            temp16alen = scale_expansion(bytcalen, bytca, bdytail, temp16a);
            temp16blen = scale_expansion(catlen, cat, bdytail, temp16b);
            temp32alen = scale_expansion(temp16blen, temp16b, 2 * bdy, temp32a);
            temp48len = expansion_sum(temp16alen, temp16a, temp32alen, temp32a, temp48);
            finlength = expansion_sum(finlength, finnow, temp48len, temp48, finother);
            finswap = finnow; finnow = finother; finother = finswap;

            temp32alen = scale_expansion(temp16blen, temp16b, bdytail, temp32a);
            temp8len = scale_expansion(cattlen, catt, bdytail, temp8);
            temp16alen = scale_expansion(temp8len, temp8, 2 * bdy, temp16a);
            temp16blen = scale_expansion(temp8len, temp8, bdytail, temp16b);
            temp32blen = expansion_sum(temp16alen, temp16a, temp16blen, temp16b, temp32b);
            temp64len = expansion_sum(temp32alen, temp32a, temp32blen, temp32b, temp64);
            finlength = expansion_sum(finlength, finnow, temp64len, temp64, finother);
            finswap = finnow; finnow = finother; finother = finswap;
        }
    }
    if (cdxtail !== 0 || cdytail !== 0) {
        if (adxtail !== 0 || adytail !== 0 || bdxtail !== 0 || bdytail !== 0) {
            $Two_Product_Sum(adxtail, bdy, adx, bdytail, u);
            n1 = -ady;
            n0 = -adytail;
            $Two_Product_Sum(bdxtail, n1, bdx, n0, v);
            abtlen = expansion_sum(4, u, 4, v, abt);

            $Cross_Product(adxtail, adytail, bdxtail, bdytail, abtt);
            abttlen = 4;
        } else {
            abt[0] = 0;
            abtlen = 1;
            abtt[0] = 0;
            abttlen = 1;
        }

        if (cdxtail !== 0) {
            temp16alen = scale_expansion(cxtablen, cxtab, cdxtail, temp16a);
            temp16blen = scale_expansion(abtlen, abt, cdxtail, temp16b);
            temp32alen = scale_expansion(temp16blen, temp16b, 2 * cdx, temp32a);
            temp48len = expansion_sum(temp16alen, temp16a, temp32alen, temp32a, temp48);
            finlength = expansion_sum(finlength, finnow, temp48len, temp48, finother);
            finswap = finnow; finnow = finother; finother = finswap;
            if (adytail !== 0) {
                temp16alen = scale_expansion_twice(4, bb, cdxtail, temp8, adytail, temp16a);
                finlength = expansion_sum(finlength, finnow, temp16alen, temp16a, finother);
                finswap = finnow; finnow = finother; finother = finswap;
            }
            if (bdytail !== 0) {
                temp16alen = scale_expansion_twice(4, aa, -cdxtail, temp8, bdytail, temp16a);
                finlength = expansion_sum(finlength, finnow, temp16alen, temp16a, finother);
                finswap = finnow; finnow = finother; finother = finswap;
            }

            temp32alen = scale_expansion(temp16blen, temp16b, cdxtail, temp32a);
            temp8len = scale_expansion(abttlen, abtt, cdxtail, temp8);
            temp16alen = scale_expansion(temp8len, temp8, 2 * cdx, temp16a);
            temp16blen = scale_expansion(temp8len, temp8, cdxtail, temp16b);
            temp32blen = expansion_sum(temp16alen, temp16a, temp16blen, temp16b, temp32b);
            temp64len = expansion_sum(temp32alen, temp32a, temp32blen, temp32b, temp64);
            finlength = expansion_sum(finlength, finnow, temp64len, temp64, finother);
            finswap = finnow; finnow = finother; finother = finswap;
        }
        if (cdytail !== 0) {
            temp16alen = scale_expansion(cytablen, cytab, cdytail, temp16a);
            temp16blen = scale_expansion(abtlen, abt, cdytail, temp16b);
            temp32alen = scale_expansion(temp16blen, temp16b, 2 * cdy, temp32a);
            temp48len = expansion_sum(temp16alen, temp16a, temp32alen, temp32a, temp48);
            finlength = expansion_sum(finlength, finnow, temp48len, temp48, finother);
            finswap = finnow; finnow = finother; finother = finswap;

            temp32alen = scale_expansion(temp16blen, temp16b, cdytail, temp32a);
            temp8len = scale_expansion(abttlen, abtt, cdytail, temp8);
            temp16alen = scale_expansion(temp8len, temp8, 2 * cdy, temp16a);
            temp16blen = scale_expansion(temp8len, temp8, cdytail, temp16b);
            temp32blen = expansion_sum(temp16alen, temp16a, temp16blen, temp16b, temp32b);
            temp64len = expansion_sum(temp32alen, temp32a, temp32blen, temp32b, temp64);
            finlength = expansion_sum(finlength, finnow, temp64len, temp64, finother);
            finswap = finnow; finnow = finother; finother = finswap;
        }
    }

    return finnow[finlength - 1];
}

export function incircle(ax, ay, bx, by, cx, cy, dx, dy) {
    const adx = ax - dx;
    const bdx = bx - dx;
    const cdx = cx - dx;
    const ady = ay - dy;
    const bdy = by - dy;
    const cdy = cy - dy;

    const bdxcdy = bdx * cdy;
    const cdxbdy = cdx * bdy;
    const alift = adx * adx + ady * ady;

    const cdxady = cdx * ady;
    const adxcdy = adx * cdy;
    const blift = bdx * bdx + bdy * bdy;

    const adxbdy = adx * bdy;
    const bdxady = bdx * ady;
    const clift = cdx * cdx + cdy * cdy;

    const det =
        alift * (bdxcdy - cdxbdy) +
        blift * (cdxady - adxcdy) +
        clift * (adxbdy - bdxady);

    const permanent =
        (Math.abs(bdxcdy) + Math.abs(cdxbdy)) * alift +
        (Math.abs(cdxady) + Math.abs(adxcdy)) * blift +
        (Math.abs(adxbdy) + Math.abs(bdxady)) * clift;

    const errbound = iccerrboundA * permanent;

    if (det > errbound || -det > errbound) {
        return det;
    }
    return incircleadapt(ax, ay, bx, by, cx, cy, dx, dy, permanent);
}

export function incirclefast(ax, ay, bx, by, cx, cy, dx, dy) {
    const adx = ax - dx;
    const ady = ay - dy;
    const bdx = bx - dx;
    const bdy = by - dy;
    const cdx = cx - dx;
    const cdy = cy - dy;

    const abdet = adx * bdy - bdx * ady;
    const bcdet = bdx * cdy - cdx * bdy;
    const cadet = cdx * ady - adx * cdy;
    const alift = adx * adx + ady * ady;
    const blift = bdx * bdx + bdy * bdy;
    const clift = cdx * cdx + cdy * cdy;

    return alift * bcdet + blift * cadet + clift * abdet;
}

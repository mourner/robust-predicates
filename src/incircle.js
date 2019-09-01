
import {
    epsilon, splitter, resulterrbound, estimate,
    fast_expansion_sum_zeroelim, scale_expansion_zeroelim
} from './util.js';

const iccerrboundA = (10.0 + 96.0 * epsilon) * epsilon;
const iccerrboundB = (4.0 + 48.0 * epsilon) * epsilon;
const iccerrboundC = (44.0 + 576.0 * epsilon) * epsilon * epsilon;

const bc = new Float64Array(4);
const ca = new Float64Array(4);
const ab = new Float64Array(4);
const axbc = new Float64Array(8);
const axxbc = new Float64Array(16);
const aybc = new Float64Array(8);
const ayybc = new Float64Array(16);
const adet = new Float64Array(32);
const bxca = new Float64Array(8);
const bxxca = new Float64Array(16);
const byca = new Float64Array(8);
const byyca = new Float64Array(16);
const bdet = new Float64Array(32);
const cxab = new Float64Array(8);
const cxxab = new Float64Array(16);
const cyab = new Float64Array(8);
const cyyab = new Float64Array(16);
const cdet = new Float64Array(32);
const abdet = new Float64Array(64);
const fin1 = new Float64Array(1152);
const fin2 = new Float64Array(1152);
const aa = new Float64Array(4);
const bb = new Float64Array(4);
const cc = new Float64Array(4);
const u = new Float64Array(4);
const v = new Float64Array(4);
const temp8 = new Float64Array(8);
const temp16a = new Float64Array(16);
const temp16b = new Float64Array(16);
const temp16c = new Float64Array(16);
const temp32a = new Float64Array(32);
const temp32b = new Float64Array(32);
const temp48 = new Float64Array(48);
const temp64 = new Float64Array(64);
const axtbb = new Float64Array(8);
const axtcc = new Float64Array(8);
const aytbb = new Float64Array(8);
const aytcc = new Float64Array(8);
const bxtaa = new Float64Array(8);
const bxtcc = new Float64Array(8);
const bytaa = new Float64Array(8);
const bytcc = new Float64Array(8);
const cxtaa = new Float64Array(8);
const cxtbb = new Float64Array(8);
const cytaa = new Float64Array(8);
const cytbb = new Float64Array(8);
const axtbc = new Float64Array(8);
const aytbc = new Float64Array(8);
const bxtca = new Float64Array(8);
const bytca = new Float64Array(8);
const cxtab = new Float64Array(8);
const cytab = new Float64Array(8);
const axtbct = new Float64Array(16);
const aytbct = new Float64Array(16);
const bxtcat = new Float64Array(16);
const bytcat = new Float64Array(16);
const cxtabt = new Float64Array(16);
const cytabt = new Float64Array(16);
const axtbctt = new Float64Array(8);
const aytbctt = new Float64Array(8);
const bxtcatt = new Float64Array(8);
const bytcatt = new Float64Array(8);
const cxtabtt = new Float64Array(8);
const cytabtt = new Float64Array(8);
const abt = new Float64Array(8);
const bct = new Float64Array(8);
const cat = new Float64Array(8);
const abtt = new Float64Array(4);
const bctt = new Float64Array(4);
const catt = new Float64Array(4);

function incircleadapt(pax, pay, pbx, pby, pcx, pcy, pdx, pdy, permanent) {
    let bdxcdy1, cdxbdy1, cdxady1, adxcdy1, adxbdy1, bdxady1;
    let bdxcdy0, cdxbdy0, cdxady0, adxcdy0, adxbdy0, bdxady0;
    let bc3, ca3, ab3;
    let finnow, finother, finswap, finlength;

    let adxtail, bdxtail, cdxtail, adytail, bdytail, cdytail;
    let adxadx1, adyady1, bdxbdx1, bdybdy1, cdxcdx1, cdycdy1;
    let adxadx0, adyady0, bdxbdx0, bdybdy0, cdxcdx0, cdycdy0;
    let aa3, bb3, cc3;
    let ti1, tj1;
    let ti0, tj0;
    let u3, v3;
    let temp8len, temp16alen, temp16blen, temp16clen;
    let temp32alen, temp32blen, temp48len, temp64len;
    let axtbblen, axtcclen, aytbblen, aytcclen;
    let bxtaalen, bxtcclen, bytaalen, bytcclen;
    let cxtaalen, cxtbblen, cytaalen, cytbblen;
    let axtbclen, aytbclen, bxtcalen, bytcalen, cxtablen, cytablen;
    let axtbctlen, aytbctlen, bxtcatlen, bytcatlen, cxtabtlen, cytabtlen;
    let axtbcttlen, aytbcttlen, bxtcattlen, bytcattlen, cxtabttlen, cytabttlen;
    let abtlen, bctlen, catlen;
    let abttlen, bcttlen, cattlen;
    let abtt3, bctt3, catt3;
    let negate;

    let bvirt, c, ahi, alo, bhi, blo, _i, _j, _0;

    const adx = pax - pdx;
    const bdx = pbx - pdx;
    const cdx = pcx - pdx;
    const ady = pay - pdy;
    const bdy = pby - pdy;
    const cdy = pcy - pdy;

    $Two_Product(bdx, cdy, bdxcdy1, bdxcdy0);
    $Two_Product(cdx, bdy, cdxbdy1, cdxbdy0);
    $Two_Two_Diff(bdxcdy1, bdxcdy0, cdxbdy1, cdxbdy0, bc3, bc[2], bc[1], bc[0]);
    bc[3] = bc3;
    const axbclen = scale_expansion_zeroelim(4, bc, adx, axbc);
    const axxbclen = scale_expansion_zeroelim(axbclen, axbc, adx, axxbc);
    const aybclen = scale_expansion_zeroelim(4, bc, ady, aybc);
    const ayybclen = scale_expansion_zeroelim(aybclen, aybc, ady, ayybc);
    const alen = fast_expansion_sum_zeroelim(axxbclen, axxbc, ayybclen, ayybc, adet);

    $Two_Product(cdx, ady, cdxady1, cdxady0);
    $Two_Product(adx, cdy, adxcdy1, adxcdy0);
    $Two_Two_Diff(cdxady1, cdxady0, adxcdy1, adxcdy0, ca3, ca[2], ca[1], ca[0]);
    ca[3] = ca3;
    const bxcalen = scale_expansion_zeroelim(4, ca, bdx, bxca);
    const bxxcalen = scale_expansion_zeroelim(bxcalen, bxca, bdx, bxxca);
    const bycalen = scale_expansion_zeroelim(4, ca, bdy, byca);
    const byycalen = scale_expansion_zeroelim(bycalen, byca, bdy, byyca);
    const blen = fast_expansion_sum_zeroelim(bxxcalen, bxxca, byycalen, byyca, bdet);

    $Two_Product(adx, bdy, adxbdy1, adxbdy0);
    $Two_Product(bdx, ady, bdxady1, bdxady0);
    $Two_Two_Diff(adxbdy1, adxbdy0, bdxady1, bdxady0, ab3, ab[2], ab[1], ab[0]);
    ab[3] = ab3;
    const cxablen = scale_expansion_zeroelim(4, ab, cdx, cxab);
    const cxxablen = scale_expansion_zeroelim(cxablen, cxab, cdx, cxxab);
    const cyablen = scale_expansion_zeroelim(4, ab, cdy, cyab);
    const cyyablen = scale_expansion_zeroelim(cyablen, cyab, cdy, cyyab);
    const clen = fast_expansion_sum_zeroelim(cxxablen, cxxab, cyyablen, cyyab, cdet);

    const ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, abdet);
    finlength = fast_expansion_sum_zeroelim(ablen, abdet, clen, cdet, fin1);

    let det = estimate(finlength, fin1);
    let errbound = iccerrboundB * permanent;
    if ((det >= errbound) || (-det >= errbound)) {
        return det;
    }

    $Two_Diff_Tail(pax, pdx, adx, adxtail);
    $Two_Diff_Tail(pay, pdy, ady, adytail);
    $Two_Diff_Tail(pbx, pdx, bdx, bdxtail);
    $Two_Diff_Tail(pby, pdy, bdy, bdytail);
    $Two_Diff_Tail(pcx, pdx, cdx, cdxtail);
    $Two_Diff_Tail(pcy, pdy, cdy, cdytail);
    if (adxtail === 0 && bdxtail === 0 && cdxtail === 0 && adytail === 0 && bdytail === 0 && cdytail === 0) {
        return det;
    }

    errbound = iccerrboundC * permanent + resulterrbound * Math.abs(det);
    det += ((adx * adx + ady * ady) * ((bdx * cdytail + cdy * bdxtail) - (bdy * cdxtail + cdx * bdytail)) +
        2.0 * (adx * adxtail + ady * adytail) * (bdx * cdy - bdy * cdx)) +
        ((bdx * bdx + bdy * bdy) * ((cdx * adytail + ady * cdxtail) - (cdy * adxtail + adx * cdytail)) +
        2.0 * (bdx * bdxtail + bdy * bdytail) * (cdx * ady - cdy * adx)) +
        ((cdx * cdx + cdy * cdy) * ((adx * bdytail + bdy * adxtail) - (ady * bdxtail + bdx * adytail)) +
        2.0 * (cdx * cdxtail + cdy * cdytail) * (adx * bdy - ady * bdx));

    if ((det >= errbound) || (-det >= errbound)) {
        return det;
    }

    finnow = fin1;
    finother = fin2;

    if ((bdxtail !== 0.0) || (bdytail !== 0.0) || (cdxtail !== 0.0) || (cdytail !== 0.0)) {
        $Square(adx, adxadx1, adxadx0);
        $Square(ady, adyady1, adyady0);
        $Two_Two_Sum(adxadx1, adxadx0, adyady1, adyady0, aa3, aa[2], aa[1], aa[0]);
        aa[3] = aa3;
    }
    if ((cdxtail !== 0.0) || (cdytail !== 0.0) || (adxtail !== 0.0) || (adytail !== 0.0)) {
        $Square(bdx, bdxbdx1, bdxbdx0);
        $Square(bdy, bdybdy1, bdybdy0);
        $Two_Two_Sum(bdxbdx1, bdxbdx0, bdybdy1, bdybdy0, bb3, bb[2], bb[1], bb[0]);
        bb[3] = bb3;
    }
    if ((adxtail !== 0.0) || (adytail !== 0.0) || (bdxtail !== 0.0) || (bdytail !== 0.0)) {
        $Square(cdx, cdxcdx1, cdxcdx0);
        $Square(cdy, cdycdy1, cdycdy0);
        $Two_Two_Sum(cdxcdx1, cdxcdx0, cdycdy1, cdycdy0, cc3, cc[2], cc[1], cc[0]);
        cc[3] = cc3;
    }

    if (adxtail !== 0.0) {
        axtbclen = scale_expansion_zeroelim(4, bc, adxtail, axtbc);
        temp16alen = scale_expansion_zeroelim(axtbclen, axtbc, 2.0 * adx, temp16a);

        axtcclen = scale_expansion_zeroelim(4, cc, adxtail, axtcc);
        temp16blen = scale_expansion_zeroelim(axtcclen, axtcc, bdy, temp16b);

        axtbblen = scale_expansion_zeroelim(4, bb, adxtail, axtbb);
        temp16clen = scale_expansion_zeroelim(axtbblen, axtbb, -cdy, temp16c);

        temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp16blen, temp16b, temp32a);
        temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c, temp32alen, temp32a, temp48);
        finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, temp48, finother);
        finswap = finnow; finnow = finother; finother = finswap;
    }
    if (adytail !== 0.0) {
        aytbclen = scale_expansion_zeroelim(4, bc, adytail, aytbc);
        temp16alen = scale_expansion_zeroelim(aytbclen, aytbc, 2.0 * ady, temp16a);

        aytbblen = scale_expansion_zeroelim(4, bb, adytail, aytbb);
        temp16blen = scale_expansion_zeroelim(aytbblen, aytbb, cdx, temp16b);

        aytcclen = scale_expansion_zeroelim(4, cc, adytail, aytcc);
        temp16clen = scale_expansion_zeroelim(aytcclen, aytcc, -bdx, temp16c);

        temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp16blen, temp16b, temp32a);
        temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c, temp32alen, temp32a, temp48);
        finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, temp48, finother);
        finswap = finnow; finnow = finother; finother = finswap;
    }
    if (bdxtail !== 0.0) {
        bxtcalen = scale_expansion_zeroelim(4, ca, bdxtail, bxtca);
        temp16alen = scale_expansion_zeroelim(bxtcalen, bxtca, 2.0 * bdx, temp16a);

        bxtaalen = scale_expansion_zeroelim(4, aa, bdxtail, bxtaa);
        temp16blen = scale_expansion_zeroelim(bxtaalen, bxtaa, cdy, temp16b);

        bxtcclen = scale_expansion_zeroelim(4, cc, bdxtail, bxtcc);
        temp16clen = scale_expansion_zeroelim(bxtcclen, bxtcc, -ady, temp16c);

        temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp16blen, temp16b, temp32a);
        temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c, temp32alen, temp32a, temp48);
        finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, temp48, finother);
        finswap = finnow; finnow = finother; finother = finswap;
    }
    if (bdytail !== 0.0) {
        bytcalen = scale_expansion_zeroelim(4, ca, bdytail, bytca);
        temp16alen = scale_expansion_zeroelim(bytcalen, bytca, 2.0 * bdy, temp16a);

        bytcclen = scale_expansion_zeroelim(4, cc, bdytail, bytcc);
        temp16blen = scale_expansion_zeroelim(bytcclen, bytcc, adx, temp16b);

        bytaalen = scale_expansion_zeroelim(4, aa, bdytail, bytaa);
        temp16clen = scale_expansion_zeroelim(bytaalen, bytaa, -cdx, temp16c);

        temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp16blen, temp16b, temp32a);
        temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c, temp32alen, temp32a, temp48);
        finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, temp48, finother);
        finswap = finnow; finnow = finother; finother = finswap;
    }
    if (cdxtail !== 0.0) {
        cxtablen = scale_expansion_zeroelim(4, ab, cdxtail, cxtab);
        temp16alen = scale_expansion_zeroelim(cxtablen, cxtab, 2.0 * cdx, temp16a);

        cxtbblen = scale_expansion_zeroelim(4, bb, cdxtail, cxtbb);
        temp16blen = scale_expansion_zeroelim(cxtbblen, cxtbb, ady, temp16b);

        cxtaalen = scale_expansion_zeroelim(4, aa, cdxtail, cxtaa);
        temp16clen = scale_expansion_zeroelim(cxtaalen, cxtaa, -bdy, temp16c);

        temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp16blen, temp16b, temp32a);
        temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c, temp32alen, temp32a, temp48);
        finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, temp48, finother);
        finswap = finnow; finnow = finother; finother = finswap;
    }
    if (cdytail !== 0.0) {
        cytablen = scale_expansion_zeroelim(4, ab, cdytail, cytab);
        temp16alen = scale_expansion_zeroelim(cytablen, cytab, 2.0 * cdy, temp16a);

        cytaalen = scale_expansion_zeroelim(4, aa, cdytail, cytaa);
        temp16blen = scale_expansion_zeroelim(cytaalen, cytaa, bdx, temp16b);

        cytbblen = scale_expansion_zeroelim(4, bb, cdytail, cytbb);
        temp16clen = scale_expansion_zeroelim(cytbblen, cytbb, -adx, temp16c);

        temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp16blen, temp16b, temp32a);
        temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c, temp32alen, temp32a, temp48);
        finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, temp48, finother);
        finswap = finnow; finnow = finother; finother = finswap;
    }

    if ((adxtail !== 0.0) || (adytail !== 0.0)) {
        if ((bdxtail !== 0.0) || (bdytail !== 0.0) || (cdxtail !== 0.0) || (cdytail !== 0.0)) {
            $Two_Product(bdxtail, cdy, ti1, ti0);
            $Two_Product(bdx, cdytail, tj1, tj0);
            $Two_Two_Sum(ti1, ti0, tj1, tj0, u3, u[2], u[1], u[0]);
            u[3] = u3;
            negate = -bdy;
            $Two_Product(cdxtail, negate, ti1, ti0);
            negate = -bdytail;
            $Two_Product(cdx, negate, tj1, tj0);
            $Two_Two_Sum(ti1, ti0, tj1, tj0, v3, v[2], v[1], v[0]);
            v[3] = v3;
            bctlen = fast_expansion_sum_zeroelim(4, u, 4, v, bct);

            $Two_Product(bdxtail, cdytail, ti1, ti0);
            $Two_Product(cdxtail, bdytail, tj1, tj0);
            $Two_Two_Diff(ti1, ti0, tj1, tj0, bctt3, bctt[2], bctt[1], bctt[0]);
            bctt[3] = bctt3;
            bcttlen = 4;
        } else {
            bct[0] = 0.0;
            bctlen = 1;
            bctt[0] = 0.0;
            bcttlen = 1;
        }

        if (adxtail !== 0.0) {
            temp16alen = scale_expansion_zeroelim(axtbclen, axtbc, adxtail, temp16a);
            axtbctlen = scale_expansion_zeroelim(bctlen, bct, adxtail, axtbct);
            temp32alen = scale_expansion_zeroelim(axtbctlen, axtbct, 2.0 * adx, temp32a);
            temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp32alen, temp32a, temp48);
            finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, temp48, finother);
            finswap = finnow; finnow = finother; finother = finswap;

            if (bdytail !== 0.0) {
                temp8len = scale_expansion_zeroelim(4, cc, adxtail, temp8);
                temp16alen = scale_expansion_zeroelim(temp8len, temp8, bdytail, temp16a);
                finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen, temp16a, finother);
                finswap = finnow; finnow = finother; finother = finswap;
            }
            if (cdytail !== 0.0) {
                temp8len = scale_expansion_zeroelim(4, bb, -adxtail, temp8);
                temp16alen = scale_expansion_zeroelim(temp8len, temp8, cdytail, temp16a);
                finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen, temp16a, finother);
                finswap = finnow; finnow = finother; finother = finswap;
            }

            temp32alen = scale_expansion_zeroelim(axtbctlen, axtbct, adxtail, temp32a);
            axtbcttlen = scale_expansion_zeroelim(bcttlen, bctt, adxtail, axtbctt);
            temp16alen = scale_expansion_zeroelim(axtbcttlen, axtbctt, 2.0 * adx, temp16a);
            temp16blen = scale_expansion_zeroelim(axtbcttlen, axtbctt, adxtail, temp16b);
            temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp16blen, temp16b, temp32b);
            temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a, temp32blen, temp32b, temp64);
            finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len, temp64, finother);
            finswap = finnow; finnow = finother; finother = finswap;
        }

        if (adytail !== 0.0) {
            temp16alen = scale_expansion_zeroelim(aytbclen, aytbc, adytail, temp16a);
            aytbctlen = scale_expansion_zeroelim(bctlen, bct, adytail, aytbct);
            temp32alen = scale_expansion_zeroelim(aytbctlen, aytbct, 2.0 * ady, temp32a);
            temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp32alen, temp32a, temp48);
            finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, temp48, finother);
            finswap = finnow; finnow = finother; finother = finswap;

            temp32alen = scale_expansion_zeroelim(aytbctlen, aytbct, adytail, temp32a);
            aytbcttlen = scale_expansion_zeroelim(bcttlen, bctt, adytail, aytbctt);
            temp16alen = scale_expansion_zeroelim(aytbcttlen, aytbctt, 2.0 * ady, temp16a);
            temp16blen = scale_expansion_zeroelim(aytbcttlen, aytbctt, adytail, temp16b);
            temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp16blen, temp16b, temp32b);
            temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a, temp32blen, temp32b, temp64);
            finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len, temp64, finother);
            finswap = finnow; finnow = finother; finother = finswap;
        }
    }
    if ((bdxtail !== 0.0) || (bdytail !== 0.0)) {
        if ((cdxtail !== 0.0) || (cdytail !== 0.0) || (adxtail !== 0.0) || (adytail !== 0.0)) {
            $Two_Product(cdxtail, ady, ti1, ti0);
            $Two_Product(cdx, adytail, tj1, tj0);
            $Two_Two_Sum(ti1, ti0, tj1, tj0, u3, u[2], u[1], u[0]);
            u[3] = u3;
            negate = -cdy;
            $Two_Product(adxtail, negate, ti1, ti0);
            negate = -cdytail;
            $Two_Product(adx, negate, tj1, tj0);
            $Two_Two_Sum(ti1, ti0, tj1, tj0, v3, v[2], v[1], v[0]);
            v[3] = v3;
            catlen = fast_expansion_sum_zeroelim(4, u, 4, v, cat);

            $Two_Product(cdxtail, adytail, ti1, ti0);
            $Two_Product(adxtail, cdytail, tj1, tj0);
            $Two_Two_Diff(ti1, ti0, tj1, tj0, catt3, catt[2], catt[1], catt[0]);
            catt[3] = catt3;
            cattlen = 4;
        } else {
            cat[0] = 0.0;
            catlen = 1;
            catt[0] = 0.0;
            cattlen = 1;
        }

        if (bdxtail !== 0.0) {
            temp16alen = scale_expansion_zeroelim(bxtcalen, bxtca, bdxtail, temp16a);
            bxtcatlen = scale_expansion_zeroelim(catlen, cat, bdxtail, bxtcat);
            temp32alen = scale_expansion_zeroelim(bxtcatlen, bxtcat, 2.0 * bdx, temp32a);
            temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp32alen, temp32a, temp48);
            finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, temp48, finother);
            finswap = finnow; finnow = finother; finother = finswap;
            if (cdytail !== 0.0) {
                temp8len = scale_expansion_zeroelim(4, aa, bdxtail, temp8);
                temp16alen = scale_expansion_zeroelim(temp8len, temp8, cdytail, temp16a);
                finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen, temp16a, finother);
                finswap = finnow; finnow = finother; finother = finswap;
            }
            if (adytail !== 0.0) {
                temp8len = scale_expansion_zeroelim(4, cc, -bdxtail, temp8);
                temp16alen = scale_expansion_zeroelim(temp8len, temp8, adytail, temp16a);
                finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen, temp16a, finother);
                finswap = finnow; finnow = finother; finother = finswap;
            }

            temp32alen = scale_expansion_zeroelim(bxtcatlen, bxtcat, bdxtail, temp32a);
            bxtcattlen = scale_expansion_zeroelim(cattlen, catt, bdxtail, bxtcatt);
            temp16alen = scale_expansion_zeroelim(bxtcattlen, bxtcatt, 2.0 * bdx, temp16a);
            temp16blen = scale_expansion_zeroelim(bxtcattlen, bxtcatt, bdxtail, temp16b);
            temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp16blen, temp16b, temp32b);
            temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a, temp32blen, temp32b, temp64);
            finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len, temp64, finother);
            finswap = finnow; finnow = finother; finother = finswap;
        }
        if (bdytail !== 0.0) {
            temp16alen = scale_expansion_zeroelim(bytcalen, bytca, bdytail, temp16a);
            bytcatlen = scale_expansion_zeroelim(catlen, cat, bdytail, bytcat);
            temp32alen = scale_expansion_zeroelim(bytcatlen, bytcat, 2.0 * bdy, temp32a);
            temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp32alen, temp32a, temp48);
            finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, temp48, finother);
            finswap = finnow; finnow = finother; finother = finswap;

            temp32alen = scale_expansion_zeroelim(bytcatlen, bytcat, bdytail, temp32a);
            bytcattlen = scale_expansion_zeroelim(cattlen, catt, bdytail, bytcatt);
            temp16alen = scale_expansion_zeroelim(bytcattlen, bytcatt, 2.0 * bdy, temp16a);
            temp16blen = scale_expansion_zeroelim(bytcattlen, bytcatt, bdytail, temp16b);
            temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp16blen, temp16b, temp32b);
            temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a, temp32blen, temp32b, temp64);
            finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len, temp64, finother);
            finswap = finnow; finnow = finother; finother = finswap;
        }
    }
    if ((cdxtail !== 0.0) || (cdytail !== 0.0)) {
        if ((adxtail !== 0.0) || (adytail !== 0.0) || (bdxtail !== 0.0) || (bdytail !== 0.0)) {
            $Two_Product(adxtail, bdy, ti1, ti0);
            $Two_Product(adx, bdytail, tj1, tj0);
            $Two_Two_Sum(ti1, ti0, tj1, tj0, u3, u[2], u[1], u[0]);
            u[3] = u3;
            negate = -ady;
            $Two_Product(bdxtail, negate, ti1, ti0);
            negate = -adytail;
            $Two_Product(bdx, negate, tj1, tj0);
            $Two_Two_Sum(ti1, ti0, tj1, tj0, v3, v[2], v[1], v[0]);
            v[3] = v3;
            abtlen = fast_expansion_sum_zeroelim(4, u, 4, v, abt);

            $Two_Product(adxtail, bdytail, ti1, ti0);
            $Two_Product(bdxtail, adytail, tj1, tj0);
            $Two_Two_Diff(ti1, ti0, tj1, tj0, abtt3, abtt[2], abtt[1], abtt[0]);
            abtt[3] = abtt3;
            abttlen = 4;
        } else {
            abt[0] = 0.0;
            abtlen = 1;
            abtt[0] = 0.0;
            abttlen = 1;
        }

        if (cdxtail !== 0.0) {
            temp16alen = scale_expansion_zeroelim(cxtablen, cxtab, cdxtail, temp16a);
            cxtabtlen = scale_expansion_zeroelim(abtlen, abt, cdxtail, cxtabt);
            temp32alen = scale_expansion_zeroelim(cxtabtlen, cxtabt, 2.0 * cdx, temp32a);
            temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp32alen, temp32a, temp48);
            finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, temp48, finother);
            finswap = finnow; finnow = finother; finother = finswap;
            if (adytail !== 0.0) {
                temp8len = scale_expansion_zeroelim(4, bb, cdxtail, temp8);
                temp16alen = scale_expansion_zeroelim(temp8len, temp8, adytail, temp16a);
                finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen, temp16a, finother);
                finswap = finnow; finnow = finother; finother = finswap;
            }
            if (bdytail !== 0.0) {
                temp8len = scale_expansion_zeroelim(4, aa, -cdxtail, temp8);
                temp16alen = scale_expansion_zeroelim(temp8len, temp8, bdytail, temp16a);
                finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen, temp16a, finother);
                finswap = finnow; finnow = finother; finother = finswap;
            }

            temp32alen = scale_expansion_zeroelim(cxtabtlen, cxtabt, cdxtail, temp32a);
            cxtabttlen = scale_expansion_zeroelim(abttlen, abtt, cdxtail, cxtabtt);
            temp16alen = scale_expansion_zeroelim(cxtabttlen, cxtabtt, 2.0 * cdx, temp16a);
            temp16blen = scale_expansion_zeroelim(cxtabttlen, cxtabtt, cdxtail, temp16b);
            temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp16blen, temp16b, temp32b);
            temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a, temp32blen, temp32b, temp64);
            finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len, temp64, finother);
            finswap = finnow; finnow = finother; finother = finswap;
        }
        if (cdytail !== 0.0) {
            temp16alen = scale_expansion_zeroelim(cytablen, cytab, cdytail, temp16a);
            cytabtlen = scale_expansion_zeroelim(abtlen, abt, cdytail, cytabt);
            temp32alen = scale_expansion_zeroelim(cytabtlen, cytabt, 2.0 * cdy, temp32a);
            temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp32alen, temp32a, temp48);
            finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, temp48, finother);
            finswap = finnow; finnow = finother; finother = finswap;

            temp32alen = scale_expansion_zeroelim(cytabtlen, cytabt, cdytail, temp32a);
            cytabttlen = scale_expansion_zeroelim(abttlen, abtt, cdytail, cytabtt);
            temp16alen = scale_expansion_zeroelim(cytabttlen, cytabtt, 2.0 * cdy, temp16a);
            temp16blen = scale_expansion_zeroelim(cytabttlen, cytabtt, cdytail, temp16b);
            temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp16blen, temp16b, temp32b);
            temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a, temp32blen, temp32b, temp64);
            finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len, temp64, finother);
            finswap = finnow; finnow = finother; finother = finswap;
        }
    }

    return finnow[finlength - 1];
}

export function incircle(pax, pay, pbx, pby, pcx, pcy, pdx, pdy) {
    const adx = pax - pdx;
    const bdx = pbx - pdx;
    const cdx = pcx - pdx;
    const ady = pay - pdy;
    const bdy = pby - pdy;
    const cdy = pcy - pdy;

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

    if ((det > errbound) || (-det > errbound)) {
        return det;
    }

    return incircleadapt(pax, pay, pbx, pby, pcx, pcy, pdx, pdy, permanent);
}

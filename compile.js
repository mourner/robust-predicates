
const fs = require('fs');
const path = require('path');

const macros = {};

macros.Fast_Two_Sum = (a, b, x, y) => `
    ${x} = ${a} + ${b};
    ${y} = ${b} - (${x} - ${a});`;

macros.Two_Sum = (a, b, x, y) => `
    ${x} = ${a} + ${b};
    bvirt = ${x} - ${a};
    ${y} = ${a} - (${x} - bvirt) + (${b} - bvirt);`;

macros.Two_Diff_Tail = (a, b, x, y) => `
    bvirt = ${a} - ${x};
    ${y} = ${a} - (${x} + bvirt) + (bvirt - ${b});`;

macros.Two_Diff = (a, b, x, y) => `
    ${x} = ${a} - ${b};
    ${macros.Two_Diff_Tail(a, b, x, y)}`;

macros.Split = (a, ahi, alo) => `
    c = splitter * ${a};
    ${ahi} = c - (c - ${a});
    ${alo} = ${a} - ${ahi};`;

macros.Two_Product = (a, b, x, y) => `
    ${x} = ${a} * ${b};
    ${macros.Split(a, 'ahi', 'alo')}
    ${macros.Split(b, 'bhi', 'blo')}
    ${y} = alo * blo - (${x} - ahi * bhi - alo * bhi - ahi * blo);`;

macros.Two_Product_Presplit = (a, b, bhi, blo, x, y) => `
    ${x} = ${a} * ${b};
    ${macros.Split(a, 'ahi', 'alo')}
    ${y} = alo * ${blo} - (${x} - ahi * ${bhi} - alo * ${bhi} - ahi * ${blo});`;

macros.Square = (a, x, y) => `
    ${x} = ${a} * ${a};
    ${macros.Split(a, 'ahi', 'alo')}
    ${y} = alo * alo - (${x} - ahi * ahi - (ahi + ahi) * alo);`;

macros.Two_One_Sum = (a1, a0, b, x2, x1, x0) => `
    ${macros.Two_Sum(a0, b, '_i', x0)}
    ${macros.Two_Sum(a1, '_i', x2, x1)}`;

macros.Two_One_Diff = (a1, a0, b, x2, x1, x0) => `
    ${macros.Two_Diff(a0, b, '_i', x0)}
    ${macros.Two_Sum(a1, '_i', x2, x1)}`;

macros.Two_Two_Sum = (a1, a0, b1, b0, x3, x2, x1, x0) => `
    ${macros.Two_One_Sum(a1, a0, b0, '_j', '_0', x0)}
    ${macros.Two_One_Sum('_j', '_0', b1, x3, x2, x1)}`;

macros.Two_Two_Diff = (a1, a0, b1, b0, x3, x2, x1, x0) => `
    ${macros.Two_One_Diff(a1, a0, b0, '_j', '_0', x0)}
    ${macros.Two_One_Diff('_j', '_0', b1, x3, x2, x1)}`;

macros.Two_One_Product = (a1, a0, b, D) => `
    ${macros.Split(b, 'bhi', 'blo')}
    ${macros.Two_Product_Presplit(a0, b, 'bhi', 'blo', '_i', `${D}[0]`)}
    ${macros.Two_Product_Presplit(a1, b, 'bhi', 'blo', '_j', '_0')}
    ${macros.Two_Sum('_i', '_0', '_k', `${D}[1]`)}
    ${macros.Fast_Two_Sum('_j', '_k', 'u3', `${D}[2]`)}
    ${D}[3] = u3;`;

macros.Cross_Product = (a, b, c, d, D, u3 = 'u3') => `
    ${macros.Two_Product(a, d, 's1', 's0')}
    ${macros.Two_Product(c, b, 't1', 't0')}
    ${macros.Two_Two_Diff('s1', 's0', 't1', 't0', u3, `${D}[2]`, `${D}[1]`, `${D}[0]`)}
    ${D}[3] = ${u3};`;

macros.Two_Product_Sum = (a, b, c, d, D) => `
    ${macros.Two_Product(a, b, 's1', 's0')}
    ${macros.Two_Product(c, d, 't1', 't0')}
    ${macros.Two_Two_Sum('s1', 's0', 't1', 't0', 'u3', `${D}[2]`, `${D}[1]`, `${D}[0]`)}
    ${D}[3] = u3;`;

macros.Square_Sum = (a, b, D) => `
    ${macros.Square(a, 's1', 's0')}
    ${macros.Square(b, 't1', 't0')}
    ${macros.Two_Two_Sum('s1', 's0', 't1', 't0', 'u3', `${D}[2]`, `${D}[1]`, `${D}[0]`)}
    ${D}[3] = u3;`;

const replaceMacros = (_, indent, id, args) => macros[id](...args.split(/, +/))
    .split('\n')
    .filter(str => str.match(/\S/))
    .map(str => indent + str.trim())
    .join('\n');

const compile = src => src.replace(/^( *)\$(\w+)\((.+)\);/mg, replaceMacros);

const files = fs.readdirSync('./src');
for (const file of files) {
    const src = fs.readFileSync(path.join('./src', file), 'utf8');
    fs.writeFileSync(path.join('./esm', file), compile(src));
}

import terser from '@rollup/plugin-terser';

const output = (input, file, plugins) => ({
    input,
    output: {
        name: 'predicates',
        format: 'umd',
        indent: false,
        file,
    },
    plugins
});

const builds = (input, name) => [
    output(input, `umd/${name}.js`, []),
    output(input, `umd/${name}.min.js`, [terser()])
];

export default [
    ...builds('index.js', 'predicates'),
    ...builds('esm/orient2d.js', 'orient2d'),
    ...builds('esm/incircle.js', 'incircle'),
    ...builds('esm/orient3d.js', 'orient3d'),
    ...builds('esm/insphere.js', 'insphere')
];

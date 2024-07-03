import config from 'eslint-config-mourner';

export default [
    ...config,
    {
        rules: {
            'camelcase': 0,
            'new-cap': 0,
            'no-unused-vars': [2, {varsIgnorePattern: 'splitter|bvirt|c|[ab]hi|[ab]lo|_[ijk0]|u3|[st][01]'}],
            'no-lonely-if': 0
        },
        languageOptions: {
            globals: {
                '$Fast_Two_Sum': false,
                '$Two_Sum': false,
                '$Two_Diff_Tail': false,
                '$Split': false,
                '$Two_Product': false,
                '$Two_Product_Presplit': false,
                '$Two_One_Product': false,
                '$Cross_Product': false,
                '$Square_Sum': false,
                '$Two_Product_Sum': false
            }
        }
    }
];

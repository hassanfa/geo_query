# Changelog

## [4.4.4](https://github.com/hassanfa/geo_query/compare/v4.4.3...v4.4.4) (2024-12-10)


### Bug Fixes

* print table if fwrite is not enabled ([fb3cf25](https://github.com/hassanfa/geo_query/commit/fb3cf25518b7266d86e707a1f50e047ced3cb409))

## [4.4.3](https://github.com/hassanfa/geo_query/compare/v4.4.2...v4.4.3) (2024-12-09)


### Bug Fixes

* handle empty queries gracefully! ([5993c9d](https://github.com/hassanfa/geo_query/commit/5993c9d56b3d2138e6330c8fb1be355dc032c63c))

## [4.4.2](https://github.com/hassanfa/geo_query/compare/v4.4.1...v4.4.2) (2024-12-08)


### Bug Fixes

* remove faulty prefix from start of GSE column ([d03f490](https://github.com/hassanfa/geo_query/commit/d03f490d71376695710cc34465d8b0416d36d563))

## [4.4.1](https://github.com/hassanfa/geo_query/compare/v4.4.0...v4.4.1) (2024-12-08)


### Bug Fixes

* nested DF for GSM values ([f1990b0](https://github.com/hassanfa/geo_query/commit/f1990b0fa32c239e0639b631f0f293168483cac8))

## [4.4.0](https://github.com/hassanfa/geo_query/compare/v4.3.0...v4.4.0) (2024-12-07)


### Features

* add FTP link to the table ([0d1909b](https://github.com/hassanfa/geo_query/commit/0d1909b77436c1af49abfe185b97f55077be4a31))

## [4.3.0](https://github.com/hassanfa/geo_query/compare/v4.2.0...v4.3.0) (2024-12-06)


### Features

* add support for any sample ([46b5852](https://github.com/hassanfa/geo_query/commit/46b5852becfd1d28da5fe32549ed71c555b5a158))

## [4.2.0](https://github.com/hassanfa/geo_query/compare/v4.1.0...v4.2.0) (2024-12-05)


### Features

* add pyarrow as dependecnies ([105dcfe](https://github.com/hassanfa/geo_query/commit/105dcfe141ef02f55d24dae4edd1b03d2e169855))
* support writing output as parquet file ([920262d](https://github.com/hassanfa/geo_query/commit/920262d9248818dbc36cbafd7567d3e5c8191701))

## [4.1.0](https://github.com/hassanfa/geo_query/compare/v4.0.0...v4.1.0) (2024-12-05)


### Features

* add funcionality to handle the output file ([1aa2ba9](https://github.com/hassanfa/geo_query/commit/1aa2ba96bdc00b075f21852ae4a4caf21df8d378))
* add support to write CSV output ([90584e6](https://github.com/hassanfa/geo_query/commit/90584e66c260d315ad0ffad41e6319ba97991356))
* support exporting the output as excel ([952e79d](https://github.com/hassanfa/geo_query/commit/952e79d11ed7452b2a5974e142daa0630322cab6))


### Bug Fixes

* move comments into write_df_.. function ([56f8b15](https://github.com/hassanfa/geo_query/commit/56f8b15f93ae7f2f4e12d5ef100b847cd396861d))
* remove redunant version file ([b7cf49f](https://github.com/hassanfa/geo_query/commit/b7cf49f1c14b2e06ff1bce96bdb915def1c32fb5))
* version import ([9cb7ea2](https://github.com/hassanfa/geo_query/commit/9cb7ea2d62d440520ca78feaf00f73c20ff9ea21))

## [4.0.0](https://github.com/hassanfa/geo_query/compare/v3.1.0...v4.0.0) (2024-12-02)


### ⚠ BREAKING CHANGES

* add file output write

### Features

* add file output write ([ac6270e](https://github.com/hassanfa/geo_query/commit/ac6270ec5a5d0159da5001bd31438ec9a1371521))
* add version to init ([b7be038](https://github.com/hassanfa/geo_query/commit/b7be0388f70d4915dd0c46ae0829e8caf9f32866))


### Bug Fixes

* remove redunant version ([20d1307](https://github.com/hassanfa/geo_query/commit/20d13079ec12bfa3d48d0368e6263b06e056701e))

## [3.1.0](https://github.com/hassanfa/geo_query/compare/v3.0.0...v3.1.0) (2024-12-02)


### Features

* encapsulate all entrez functionality in a single class ([26d602a](https://github.com/hassanfa/geo_query/commit/26d602adae2c9c3a1fd784b1af3da17923683d72))

## [3.0.0](https://github.com/hassanfa/geo_query/compare/v2.1.0...v3.0.0) (2024-12-02)


### ⚠ BREAKING CHANGES

* limit output only for GSE and GSM and print df for both
* add EntrezGDS class and print an output in df format

### Features

* add EntrezGDS class and print an output in df format ([7d328b9](https://github.com/hassanfa/geo_query/commit/7d328b9bd0118f8a6337da5500d43cc36af19d46))
* add option for multiple organism ([c93e73c](https://github.com/hassanfa/geo_query/commit/c93e73c70836b06d984075201481b4bb0ff58784))
* add option to print ids of result records ([3c3c1b5](https://github.com/hassanfa/geo_query/commit/3c3c1b58a4ccd0b0927f7e21924c3051106da232))
* add polars to the requirements ([e49a353](https://github.com/hassanfa/geo_query/commit/e49a353d3ed52c6117d2899b4d171e58b6882cb6))
* limit output only for GSE and GSM and print df for both ([cde22e0](https://github.com/hassanfa/geo_query/commit/cde22e0848d41d7a8f54eaee3be8a4261aba49aa))
* move geo esearch into its own function ([2b2e1e0](https://github.com/hassanfa/geo_query/commit/2b2e1e0f5a7f18133baca009b488d56b4e1b74d8))

## [2.1.0](https://github.com/hassanfa/geo_query/compare/v2.0.0...v2.1.0) (2024-11-28)


### Features

* add count and mesh operator functionality ([79b43bb](https://github.com/hassanfa/geo_query/commit/79b43bbd6dc7568b96efff0bd731a3a6472c8349))
* add logger to print output ([982fada](https://github.com/hassanfa/geo_query/commit/982fada74f9a036b10a5985e8bbf8ceb90c77bca))

## [2.0.0](https://github.com/hassanfa/geo_query/compare/v1.1.0...v2.0.0) (2024-11-27)


### ⚠ BREAKING CHANGES

* add multiple items to sample entry and add a function for building query

### Features

* add multiple items to sample entry and add a function for building query ([9ea9581](https://github.com/hassanfa/geo_query/commit/9ea95812949827eca522448da18b3bcec8ae8fb4))

## [1.1.0](https://github.com/hassanfa/geo_query/compare/v1.0.0...v1.1.0) (2024-11-26)


### Features

* options and pip installable ([d1e6be9](https://github.com/hassanfa/geo_query/commit/d1e6be99714abc1e62eea1b61633f1649a9241c7))
* simplify and add options to query ([935058f](https://github.com/hassanfa/geo_query/commit/935058ffa08dbc1551cfba4d1fda0f0a8c406e66))
* update docs ([ead580a](https://github.com/hassanfa/geo_query/commit/ead580a2cb315b6e12a114158f7e37272b478ded))

## 1.0.0 (2024-11-26)


### Features

* action and initial commit ([51a74ee](https://github.com/hassanfa/geo_query/commit/51a74ee3b14bcbbdd7290449d0f049b4e9f35fca))
* code owners ([6d2b46f](https://github.com/hassanfa/geo_query/commit/6d2b46f4fcdcecb5038e98a96c8ca6cca1bc3606))

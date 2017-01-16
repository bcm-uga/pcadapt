context("GET_OUTPUT_NAME")

test_that("output name is correct", {
  test.lfmm <- "/R.Users/luuk/3.3/pcadapt/ext.data/geno.3.pops.lfmm"
  res.lfmm <- get.output.name(test.lfmm)
  expect_match(res.lfmm, "/R.Users/luuk/3.3/pcadapt/ext.data/geno.3.pops.pcadapt") 
  
  test.vcf <- "/R.Users/luuk/3.3/pcadapt/ext.data/geno.3.pops.vcf"
  res.vcf <- get.output.name(test.vcf)
  expect_match(res.vcf, "/R.Users/luuk/3.3/pcadapt/ext.data/geno.3.pops.pcadapt")
  
  test.ped <- "/R.Users/luuk/3.3/pcadapt/ext.data/geno.3.pops.ped"
  res.ped <- get.output.name(test.ped)
  expect_match(res.ped, "/R.Users/luuk/3.3/pcadapt/ext.data/geno.3.pops.pcadapt")
  
  test.noext <- "/R.Users/ext.data/geno.3.pops"
  res.noext <- get.output.name(test.noext)
  expect_match(res.noext, "/R.Users/ext.data/geno.3.pops.pcadapt")
  
  test.nopt <- "/Users/extdata/geno3pops"
  res.nopt <- get.output.name(test.nopt)
  expect_match(res.nopt, "/Users/extdata/geno3pops.pcadapt")
})

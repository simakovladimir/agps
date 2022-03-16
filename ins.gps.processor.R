# Copyright 2015 Vladimir Simakov
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to
# deal in the Software without restriction, including without limitation the
# rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
# sell copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
# IN THE SOFTWARE.

igp.preprocess<-function(gps,ins)
{
# check arguments and variable "data" for validity
  data<-lapply(mget(names(formals())),read.table)
  data<-
    mapply(
      cbind,
      lapply(
        lapply(
          lapply(
            lapply(
              data,`[`,,1L:6L),
            apply,
            1L,
            as.function(
              alist(
                x=,
                as.numeric(
                  as.POSIXct(
                    strptime(
                      paste(x,collapse=" "),
                      "%Y %m %d %H %M %OS",
                      "GMT"),
                    "GMT"))))),
          `-`,
          min(
            sapply(
              lapply(
                data,`[`,1L,1L:6L),
              as.function(
                alist(
                  x=,
                  as.numeric(
                    as.POSIXct(
                      strptime(
                        paste(x,collapse=" "),
                        "%Y %m %d %H %M %OS",
                        "GMT"),
                      "GMT"))))))),
        matrix,
        dimnames=list(NULL,"t")),
      lapply(
        data,`[`,,-(1L:6L)),
      SIMPLIFY=FALSE,
      USE.NAMES=TRUE)
  return(data)
}
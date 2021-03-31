###############################################################################
#
# Author contact:
# zhaofei920810@gmail.com
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
###############################################################################

if (!requireNamespace("jsonlite", quietly = TRUE)) install.packages("jsonlite", update = F, ask = F)

args <- commandArgs(1)


print(args[1])
print(args[2])

dat <- jsonlite::fromJSON(txt = args[1])
info <- do.call(rbind.data.frame, dat["hits"])[3:7]
detail_list <- lapply(dat[["hits"]][["hit_hsps"]], function(x) x[1, ])
detail <- do.call(rbind.data.frame, detail_list)[-c(1, 9, 10, 12, 17, 18, 19)]
suppressWarnings(
  output <- cbind(unlist(dat["query_def"]), info, detail)
)
names(output)[1] <- "seqid"
output <- output[, c(1, 2, 3, 5, 6, 8, 10, 11, 15, 16, 17, 18)]
write.table(output, file = args[2], sep = "\t", quote = F, row.names = F)

print("done")
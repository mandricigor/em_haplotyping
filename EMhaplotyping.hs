
import Data.List(nub)
import qualified Data.Map as Map
import Control.Monad(replicateM)


haplotypes :: [Int] -> [[Int]]
haplotypes [] = []
haplotypes gen = haplotypes' gen []


--find all haplotypes of a given genotype
haplotypes' :: [Int] -> [[Int]] -> [[Int]]
haplotypes' gen bag
	| null gen = bag
	| null bag && not (null gen) = haplotypes' (init gen) (hapN $ last gen)
	| otherwise = haplotypes' (init gen) [x ++ y | x <- hapN $ last gen, y <- bag]
	
haplotypes'' x
	| countTwos x == 0 = [(x, x)]
	| otherwise = map unzip [x `pairHap` y | y <- genbin $ countTwos x]

pairHap [] _ = []
pairHap (x:hs) [] = (x, x): pairHap hs []
pairHap (x:hs) (s:str)
	| x == 1 = (1,1): pairHap hs (s:str)
	| x == 0 = (0,0): pairHap hs (s:str)
	| x == 2 && s == '0' = (0,1): pairHap hs str
	| x == 2 && s == '1' = (1,0): pairHap hs str


genbin n = take (2 ^ (n - 1)) $ replicateM n "01"


countTwos :: [Int] -> Int
countTwos [] = 0
countTwos (x:hs)
	| x == 2 = 1 + countTwos hs
	| otherwise = countTwos hs


--find all haplotypes present in the set of genotypes
allHaplotypes [] = []
allHaplotypes genotypes = nub $ foldl1 (++) (map haplotypes genotypes)


hapN :: Int -> [[Int]]
hapN 1 = [[1]]
hapN 0 = [[0]]
hapN 2 = [[1], [0]]


hapNumOrder hap [] = 0
hapNumOrder hap haps = hapNumOrder' hap haps 0

hapNumOrder' hap (h:haps) order 
	| hap /= h = hapNumOrder' hap haps (order + 1)
	| otherwise = order


hapGenMapping genHaps haplotypes = [hapGenMapping' genHap haplotypes | genHap <- genHaps]


hapGenMapping' genHap haplotypes = [hapGenMapping'' hapPair haplotypes | hapPair <- genHap]


hapGenMapping'' hapPair haplotypes = hapNumOrder (fst hapPair) haplotypes : [hapNumOrder (snd hapPair) haplotypes]



e hapGenMaps frequencies = [e' hapGenMap frequencies | hapGenMap <- hapGenMaps]

e' hapGenMap frequencies = [(hg, prod frequencies hg / sum (map (prod frequencies) hapGenMap)) | hg <- hapGenMap]


prod frequencies hg = (frequencies !! head hg) * (frequencies !! (hg !! 1))



m :: [[Int]] -> [[([Int], Double)]] -> [Double]
m gen genHapFreq = [m' k (concat genHapFreq) / (2 * fromIntegral (length gen)) | k <- [0..length (allHaplotypes gen) - 1]]



m' _ [] = 0
m' k (ghf:genHapFreqFl)
	| fst ghf == [k, k] = 2 * snd ghf + m' k genHapFreqFl
	| k `elem` fst ghf = snd ghf + m' k genHapFreqFl
	| otherwise = m' k genHapFreqFl


em'' :: [[Int]] -> [Double] -> [Double]
em'' x frequencies0 = m x (e (hapGenMapping (map haplotypes'' x) (allHaplotypes x)) frequencies0)



em' _ 0 frequencies = frequencies
em' x iterations frequencies0 = funcPow iterations (em'' x) frequencies0


funcPow 0 func = id
funcPow times func = func . funcPow (times - 1) func

em [] _ = []
em x iterations = zip (em' x iterations [1 / fromIntegral y | i <- [1..y]]) (allHaplotypes x)
	where y = length (allHaplotypes x)




# The Complementary and Integrative Health Lexicion(CIHLex)
The CIHLex is a colletion of psychological and physical complementary and integrative health(CIH) approach. The CIHLex was built by combining two CIH term list from our previous work:
1. CIH terms extracted from literature review(LRCIH)
2. Natural Medicine CIH terms (NMCIH)


## File structure
The CIH file contains the following keys:
- ```therapy```: The CIH therapy
- ```unique_concept```: the "unique concept" of the CIH therapy
- ```unique_terms```: list of unique terms that can represent the CIH therapy
- ```CIH_class```: the class of the CIH therapy, where 1=Physical, 2=Psychological, 3=Other
- ```num_UT```: the number of unique terms




## Example of data
"Alexander technique": 
{"unique_concept": "Alexander technique", 
"unique_terms": ["Alexander technique", "Alexander Proprioception", "AT", "FM Alexander Technique", "Technique Alexander"], 
"CIH_class": "1",
"num_UT": 5}

## Using the dataset
Use of this dataset is subject of the terms of the [LICENSE](CIH/LICENSE.md).

## About this Page
This is a github for a paper, 'Optimization of Green Space Allocation for Heat Reduction Considering Vulnerability of Elderly' (about to be published) <br>
In this github, there are codes and data used in the paper. The code is executed in MATLAB.

Feel free to contact me:) <br>
email : serendipity5786@snu.ac.kr

<br>

## Description of Code and Data : English version

### [Data]  
These are the 4 input datasets used in GA. They are in ASCII format with a resolution of 100m*100m.  
- allpop_100.txt: Total population data  
- fixed_100.txt: Development-restricted area data  
- lulc_100.txt: Land cover data  
- vulclass_100.txt: Population aged 65 and over data  

### [Code]  
- main.m: The primary GA code used initially  
- minmax.m: Code for calculating the min and max of each fitness function to find the global optimum

### [Function]  
- cluster.m: A function for clustering the aggregated greenery  
- extent.m: A function for calculating the area of each cluster  
- fitness_allpop.m: Fitness function for maximizing cooled POP  
- fitness_temp.m: Fitness function for maximizing cooling  
- fitness_vulclass.m: Fitness function for maximizing cooled P65+  
- initialize.m: A function for creating the initial 50 parent solutions  
- isEdgeCell.m: A function to determine if a cell belongs to the edge of a cluster  
- Normal.m: A function for normalization  
- PlotCosts.m: A function to visualize the changes in values as a graph  
- tournament_selection.m: Tournament selection function  
- whole.m: Code for identifying the indirect cooling extent

<br>

## Description of Code and Data : Korean version
### [Data]
GA에 input 된 4개의 데이터입니다. ASCII 형태의 데이터며, 100m*100m의 해상도를 가지고 있습니다.
- allpop_100.txt : 전체 인구 데이터
- fixed_100.txt : 개발 금지지역 데이터
- lulc_100.txt : 토지피복 데이터
- vulclass_100.txt : 65세 이상 인구 데이터

### [code]
- main.m : 처음 사용하는 GA 코드
- minmax : 각 fitness function의 min과 max를 구하는 코드, global optimum 구하기

### [function]
- cluster.m : 모여진 greenery들을 클러스터링 하는 함수
- extent.m : 각 클러스터의 area를 구하는 함수
- fitness_allpop.m : Maximize of cooled POP에 대한 fitness 함수
- fitness_temp.m : Maximize of cooling에 대한 fitness 함수
- fitness_vulclass.m : Maximize of cooled P65+에 대한 fitness 함수
- initialize.m : 초기 50개의 부모해를 만드는 함수
- isEdgeCell.m : 클러스터의 edge에 속하는지 판단하는함수
- Normal.m : 정규화하는 함수
- PlotCosts.m : 값의 변화를 그래프로 구현시키는 함수
- tournament_selection.m : tournament selection 함수
- whole.m : indirect cooling extent를 규명하는 코드

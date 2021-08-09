function [ resultPathDirectory] = Set_resultPathToDirectory (resultPath,cutoff,repetition)

      resultPathDirectory = strcat(resultPath,'1e',num2str(cutoff),'/rep',num2str(repetition),'/');  

end
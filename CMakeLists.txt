CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

SET(CMAKE_CXX_STANDARD 11)
#SET(CMAKE_BUILD_TYPE Debug)
#SET(CMAKE_CXX_FLAGS "-O0")
#SET(CMAKE_C_FLAGS "-O0")

PROJECT(main)

ADD_EXECUTABLE(main
    src/aux.cpp
    src/ReWeight.cpp
    src/MinimizeFCNEpsilon31.cpp
    src/MinimizeFCNSummedEnergy.cpp
    src/StaticsGroup.cpp
    src/read_data.cpp
    #src/Analysis.cpp
    #src/Analysis_ReadData.cpp
    #src/Analysis_CanvasFunctions.cpp
    #src/Analysis_EventLoopFunctions.cpp
    #src/Analysis_SensitivityFunctions.cpp
    #src/Analysis_SensitivityChisquare1.cpp
    #src/Analysis_SensitivityChisquare2.cpp
    #src/Analysis_SensitivityLL1.cpp
    #src/Analysis_SensitivityLL2.cpp
    #src/SubAnalysis.cpp
    #src/SubAnalysis_EventLoopFunctions.cpp
    #src/SubAnalysis_SensitivityFunctions.cpp
    #src/SubAnalysis_SensitivityChisquare1.cpp
    #src/SubAnalysis_SensitivityChisquare2.cpp
    #src/SubAnalysis_SensitivityLL1.cpp
    #src/SubAnalysis_SensitivityLL2.cpp
    src/main.cpp
)
ADD_EXECUTABLE(filestripper src/filestripper.cpp)
#ADD_EXECUTABLE(createtree
#    src/ReWeight.cpp
#    src/createtree.cpp
#)
#ADD_EXECUTABLE(sensitivity
#    src/sensitivity.cpp
#)
#ADD_EXECUTABLE(canvas_test src/canvas_test.cpp)

FIND_PACKAGE(ROOT REQUIRED COMPONENTS Minuit2)
MESSAGE(STATUS "CMAKE_PREFIX_PATH=" ${CMAKE_PREFIX_PATH})
MESSAGE(STATUS "ROOT_INCLUDE_DIRS=" ${ROOT_INCLUDE_DIRS})
MESSAGE(STATUS "ROOT_LIBRARIES=" ${ROOT_LIBRARIES})
LIST(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

# add the vector 3 library folder
INCLUDE_DIRECTORIES("include")
#INCLUDE_DIRECTORIES("~/Vector3")
#INCLUDE_DIRECTORIES("~/CERN-ROOT-Canvas-Factory")
#INCLUDE_DIRECTORIES("~/Numerical-String-Format")
#INCLUDE_DIRECTORIES("~/CPP-Debug")
INCLUDE_DIRECTORIES("~/Program-Arguments")


INCLUDE_DIRECTORIES(${ROOT_INCLUDE_DIRS})

TARGET_LINK_LIBRARIES(main ${ROOT_LIBRARIES})
TARGET_LINK_LIBRARIES(filestripper ${ROOT_LIBRARIES})
#TARGET_LINK_LIBRARIES(createtree ${ROOT_LIBRARIES})
#TARGET_LINK_LIBRARIES(sensitivity ${ROOT_LIBRARIES})
#TARGET_LINK_LIBRARIES(canvas_test ${ROOT_LIBRARIES})

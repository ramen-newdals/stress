cmake_minimum_required(VERSION 3.0.0 FATAL_ERROR)
project(stress VERSION 0.1.0 LANGUAGES CXX C)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

find_package( VTK REQUIRED NO_MODULE)
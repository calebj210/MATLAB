function [data] = un1997();
% function [data] = un1997();
% this functions returns a data structure containing economic and
% demographic indicators for 25 countries, taken from the UN 1997
% report.
%
% to access the structure: type for example
% >> data.country

data.country = [
    "Albania"
    "Argentina"
    "Australia"
    "Austria"
    "Benin"
    "Bolivia"
    "Brazil"
    "Camobida"
    "China"
    "Colombia"
    "Croatia"
    "El Salvador"
    "France"
    "Greece"
    "Guatemala"
    "Iran"
    "Italy"
    "Malawi"
    "Netherlands"
    "Pakistan"
    "Papua New Guinea"
    "Peru"
    "Romania"
    "USA"
    "Zimbabwe"];

data.increase = [
    1.2
    1.2
    1.1
    1
    3.2
    2.4
    1.5
    2.8
    1.1
    1.7
    -1.5
    2.2
    0.4
    0.6
    2.9
    2.3
    -0.2
    3.3
    0.7
    0.1
    1.9
    1.7
    -0.5
    1.1
    4.4];

data.life = [
    69.2
    68.6
    74.7
    73
    45.9
    57.7
    64
    50.1
    66.7
    66.4
    67.1
    63.9
    73
    75
    62.4
    67
    74.2
    45
    74.4
    60.6
    55.2
    64.1
    66.6
    72.5
    52.4];

data.imr = [
    30
    24
    7
    7
    86
    75
    58
    116
    44
    37
    9
    46
    7
    10
    48
    36
    8
    143
    7
    91
    68
    64
    23
    9
    67];

data.tfr = [
    2.9
    2.8
    1.9
    1.5
    7.1
    4.8
    2.9
    5.3
    2
    2.7
    1.7
    4
    1.7
    1.4
    5.4
    5
    1.3
    7.2
    1.6
    6.2
    5.1
    3.4
    1.5
    2.1
    5];

data.gdp = [
    659.91
    4343.04
    17529.98
    20561.88
    398.21
    812.19
    3219.22
    97.39
    341.31
    1246.87
    5400.66
    988.58
    21076.77
    6501.23
    831.81
    9129.34
    19204.92
    229.01
    18961.9
    358.59
    839.03
    1674.15
    1647.97
    21965.08
    686.75];

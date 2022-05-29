+++
title = "Parallelization"
weight = 20
+++

# (and other common options)
___
Throughout Captus' commands we provide common options that allow you limit the computer's resources available to `Captus`, change the way of running parallel tasks, and control the amount of text shown during a run:

### **`--ram`**
With this option you can specify the maximum RAM in GB that `Captus` is allowed to use. For example, if your system has 64 GB of RAM and you want to limit the use to 32.5 GB you would set the argument `--ram 32.5`.

This argument is optional, the default is **auto** (= 99% of available RAM).
___
### **`--threads`**
Similarly, you can specify the maximum number of CPU cores that `Captus` is allowed to use. If your system has 16 CPU cores but you need some cores for other analyses you could reduce it to, for example, 8 cores with `--threads 8`.

This argument is optional, the default is **auto** (= all CPU cores).
___
### **`--concurrent`**
This option sets the maximum number of tasks to run in parallel at any moment. For example, let's imagine you set `--threads 16` and `--concurrent 2`, this means that `Captus` will run only 2 tasks in parallel but each of those tasks can use up to 8 CPU cores. 

This argument is optional, the default is **auto** (the automatic adjustment varies between analysis steps).
___
### **`--debug`**
This flag enables the debugging mode, this _**disables parallelization**_ so errors can be logged to screen. If you are seeing some samples failing steps or some other unexpected behavior you can enable `--debug` and submit the error shown to the `Issues` (https://github.com/edgardomortiz/Captus/issues) section in our GitHub repository.
___
### **`--show_less`**
This flags produces less verbose screen printout. Essentially, information about each sample will not be shown (but still logged) during the run.

___
Created by [Edgardo M. Ortiz]({{< ref "../../more/credits/#edgardo-m-ortiz">}}) (06.08.2021)  
Last modified by [Edgardo M. Ortiz]({{< ref "../../more/credits/#edgardo-m-ortiz">}}) (29.05.2022)
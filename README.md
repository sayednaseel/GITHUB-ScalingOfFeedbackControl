# Scaling of Feedback Control

## Description
This repository contains the codes required to simulate the results from a research study on the scaling of perturbation response times across the size range of terrestrial mammals [publication](https://www.biorxiv.org/content/10.1101/2024.09.23.614404v1). We used computational models to simulate fast perturbation responses under feedback (PD) control and feedforward (bang-bang) control. We studied how control is affected by two limitations: muscle force capacity limits and sensorimotor time delays. We simulated two perturbation response scenarios: a swing leg repositing task after a trip (Swing task), and a posture recovery after a push (Posture task). We developed two types of models: scaled models that are parametrized with inertial, muscular and neural features to represent the size range of terrestrial mammals, and a normalized feedback control model to generate general predictions. 

<!---
![picture](FBandFFblock.jpg)
-->

The codes are distributed across the following folders:
- Swing Task-Scaled-Feedback
- Swing Task-Scaled-Feedforward
- Swing Task-Norm-Feedback
- Posture Task-Scaled-Feedback
- Posture Task-Scaled-Feedforward
- Posture Task-Norm-Feedback


## Instructions

### **Swing Task-Scaled-Feedback:** 
- **Master_SwingTask_Scaled.m**: Script to pass initial guess for controller gains to ddeSwingTask_Scaled.m, collect optimized results, perform linear fit, and graph results
    - Set run_opt=0/1 to run initial guess/run optimization
    - Set plotfig=0/1 to switch graphs off or on
- **ddeSwingTask_Scaled.m**: This function accepts initial guesses for controller gains, optimizes controller gains, and outputs simulation results for a single mass
- **Data_SwingTask.mat**: Dataset with results for all masses
- **Data_Inertialdelay_SwingTask.mat**: Dataset with inertial delay times for swing task

### **Swing Task-Scaled-Feedforward:** 
- **Master_SwingTask_ScaledFF.m**: Script to pass initial guess for tswitch to ddeSwingTask_ScaledFF.m, collect optimized results, perform linear fit, and graph results
    - Set run_opt=0/1 to run initial guess/run optimization
    - Set parms.plotfig=0/1 to switch graphs off or on
- **odeSwingTask_ScaledFF.m**: This function accepts initial guesses for Tswitch, optimizes tswitch, and outputs simulation results for a single mass
- **Data_SwingTaskFF.mat**: Dataset with results for all masses

### **Swing Task-Norm-Feedback:** 
- **SwingTaskNorm.slx**: Simulink model of normalized swing task with time delays and actuator saturation limits
- **init_Norm.m**: Set solver tolerances and error message settings
- **Landscape_SwingTask_Norm.m**: Code to simulate the normalized swing task model for initial guess controller gains. And to run a brute force search through Kp and Kd to determine settling times and % overshoot
- **Master_SwingTask_Norm.m**: Code to evaluate the relationship between normalized torque limits(tau_iso) and response time (tresp) for the normalized swing task model. For each tau_iso, run_SwingTask_Norm optimizes controller gains to find the fastest settling time without overshoot. Tries fitting different curves to the tau_iso vs tresp relationship
- **run_SwingTask_Norm.m**:Code to simulate normalized swing task model for a single tau_iso, and return settling times, overshoot, angle, angular velocity and torque profiles. Optimizes controller gains Kp and Kd to minimize settling time with 0 overshoot
- **optGains_SwingTask_Norm.mat**:Dataset with optimal controller gains that produced the fastest normalized settling times (7.09)
- **Landscape_SwingTask_Norm.mat**: Dataset with brute force search results on controller gains in normalized swing task. To view the settling time and overshoot landscapes from a brute force search of Kp and Kd (supplementary material figure S4), set plotfig=1, load this dataset to the workspace, and run the codeblock titled “plotting settling time and overshoot brute force search landscapes” in Landscape_SwingTask_Norm.m
- **TisofitST_v13.mat**:Dataset with optimal controller gains and response times for a range of tau_iso. Normalized swing task model

### **Posture Task-Scaled-Feedback**
- **Master_PostureTask_Scaled.m**: Script to pass initial guess for controller gains to ddePostureTask_Scaled.m, collect optimized results, perform linear fit, and graph results 
    - Set run_opt=0/1 to run initial guess/run optimization
    - Set plotfig=0/1 to switch graphs off or on
- **ddePostureTask_Scaled.m**: This function accepts initial guesses for controller gains, optimizes controller gains, and outputs simulation results for a single mass
- **Data_PostureTask.mat**: Dataset with results for all masses
- **Data_Inertialdelay_PostureTask.mat**: Dataset with inertial delay times for posture task

## References
- Effects of sensorimotor delays and muscle force capacity limits on the performance of feedforward and feedback control in animals of different sizes
Sayed Naseel Mohamed Thangal, Heather L. More, C. David Remy, J. Maxwell Donelan
bioRxiv 2024.09.23.614404; doi: https://doi.org/10.1101/2024.09.23.614404 
- More HL, Donelan JM. Scaling of sensorimotor delays in terrestrial mammals. Proceedings of the Royal Society B: Biological Sciences. 2018;285: 20180613. Available: http://dx.doi.org/10.1098/rspb.2018.0613
- Thangal, S.N.M., Donelan, J.M., 2020. Scaling of inertial delays in terrestrial mammals. PLoS One 15, e0217188. https://doi.org/10.1371/journal.pone.0217188
- Kilbourne, B.M., Hoffman, L.C., 2013. Scale effects between body size and limb design in quadrupedal mammals. PLoS One 8, e78392.
- Alexander, R.M., Jayes, A.S., Maloiy, G.M.O., Wathuta, E.M., 1981. Allometry of the leg muscles of mammals. J. Zool. 194, 539–552. https://doi.org/10.1111/j.1469-7998.1981.tb04600.x

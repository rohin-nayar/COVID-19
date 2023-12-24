# COVID-19 Virus Spread Simulation

## Overview

In response to the challenges posed by the COVID-19 pandemic, this project aims to simulate the spread of a virus within a defined rectangular space. The simulation involves a population of individuals moving freely, with one initially infected individual. The goal is to visualize and understand the dynamics of virus transmission within this confined space.

## Simulation Details

### Initial Setup

- The simulation takes place in a rectangular space defined by width and height parameters.
- Individuals are randomly placed within this space with initial positions specified by x and y coordinates.
- Each individual has a constant speed (V) and direction (θ) randomly generated.

### Health Status Visualization

- Healthy individuals are represented by the color green.
- Infected individuals are represented by the color red.
- Individuals in the incubation period are represented by the color orange.
- Recovered individuals are represented by the color blue.

### Virus Transmission

- Close proximity (within 2 meters) to an infected individual has a 50% chance of transmission.
- During the incubation period, close encounters with infected individuals have a 30% chance of transmission.
- The incubation period lasts for 2 days, and symptoms develop after that.

### Movement

- Individuals move based on their speed and direction.
- Velocity components in the x- and y-directions are calculated as follows:
  - u = V cos(θ)
  - v = V sin(θ)
- Position updates are calculated at each simulation step:
  - x_new = x_old + u * Δt
  - y_new = y_old + v * Δt

### Optional Challenges

- (Optional) Consider slowing down the speed for sick individuals.
- (Optional) Change direction slightly after close encounters with other individuals.

## Usage

1. Clone the repository:

   ```bash
   git clone https://github.com/yourusername/virus-spread-simulation.git
   cd virus-spread-simulation

2. Open MATLAB and navigate to the project directory.

3. Run the simulation script:

4. View the simulation results in MATLAB and observe the spread of the virus within the graphical representation of the rectangular space.

## Acknowledgments

This simulation project draws inspiration from the Washington Post's insightful simulations that visualized the spread of the coronavirus. The intent is to provide a simplified yet informative representation of virus spread dynamics within a confined space. Special thanks to the global scientific community and healthcare workers for their tireless efforts in addressing the challenges posed by the COVID-19 pandemic.

The project also recognizes the collaborative spirit of the open-source community and the valuable contributions made by developers and researchers in understanding and mitigating the impact of infectious diseases.

Contributors to this simulation project, whether direct or indirect, play a role in fostering awareness and understanding of virus transmission dynamics. Your collective efforts contribute to the ongoing dialogue surrounding public health and the challenges we face in unprecedented times.

# BiasSimulations
Simulations for demonstrating the directionality of bias from confounding, selection bias, and information bias.

In this script, we calculate the bias in the risk difference and the risk ratio for various data generating mechanisms representing uncontrolling confounding, improper control for a collider (selection bias) and missclassification error (information bias). We simulated binary values for $A$, $Y$, and other relevant nodes according to a data generating mechanism representing confounding, selection/collider bias, or information bias. For collider bias, we simulated scenarios in which there were no nodes between $A$ and $C$ and $Y$ and $C$, and in which we included additional, intermediate nodes between $A$ and $C$ and $Y$ and $C$, resulting in a greater number of signed edges. For each simulation, we generated either 1,000 or 10,000 observations (representing sample size of a study), varying the strength of the relationship between each node in the DAG, and we estimated the bias between the true risk difference and risk ratio and the risk difference and risk ratio that would be observed in the presence of bias. We repeated simulations 1,000 times per sample size and association strength, computed the average bias across the 1,000 simulations, and generated boxplots to represent the distribution of bias from each of the 1,000 simulations 

The following simulations are encoded:

## 1.0 Confounding

All examples in this section pertain to the following DAG structure: $A \rightarrow Z \leftarrow Y$, where $Z$ is a confounder.

**Scenaro 1:** A-Z and Y-Z edges are both positive

**Scenaro 2:** A-Z edge is positive and Y-Z edge is negative

**Scenaro 3:** A-Z and Y-Z edges are both negative

## 2.0 Selection bias (simple)

All examples in this section pertain to the following DAG structure: $A \leftarrow C \leftarrow Y$, where $C$ is a collider

**Scenaro 1:** A-C and Y-C edges are both positive

**Scenaro 2:** A-C edge is positive and Y-C edgs is negative

**Scenaro 3:** A-C and Y-C edges are both negative

### 2.1 Selection bias (intervening nodes)

All examples in this section pertain to the following DAG structure: $A \leftarrow U_A \rightarrow C \leftarrow U_Y \ rightarrow Y$, where $C$ is a collider, and $U_A$ and $U_Y$ are unmeasured intermediate nodes between $A$ and $C$ and $A$ and $Y$, respectively.

**Scenaro 1:** All edges are positive

**Scenaro 2:** One edge is negative, all other edges are positive.

## 3.0 Information bias (differential misclassification of the exposure)

All examples in this section pertain to the following DAG structure: $A \rightarrow A^* \leftarrow U \rightarrow Y$, where $A^*$ is the measured version of $A$ and $U$ is an unmeasured intermediate node between $A^*$ and $Y$.

**Scenaro 1:** $U$ is positively associated with $A^*$ and $Y$

**Scenaro 2:** $U$ is positively associated with $A^*$ and negatively associated with $Y$


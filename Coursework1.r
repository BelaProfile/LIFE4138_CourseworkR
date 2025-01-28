# R Coursework part 1

names <- c("Ash","Bash","Cash","Dash","Flash")  #column 1 with names which are strings
student_id <- c(20,30,40,50,60)  #column 2 with student_id which is an integer vector

dataframe <- data.frame(names, student_id)

dataframe

my_vector <- seq(from = 1, to = 1000, by = 1)
my_matrix <- matrix(1:100, nrow = 10, ncol = 10)
my_dataframe <- data.frame(name = rep(c("a", "b"), each = 10), value = 21:40)
my_list <- list(my_vector, my_matrix, my_dataframe)

my_list



#for loop
for (num in 200:300){
  print(num)
}

#using sapply 
sapply (200:300, print)

#basic method
x<- c(200:300)
print(x)


animals <- c("tiger", "lion", "badger", "fox", "rabbit", "fish", "dog", "octopus")

animals   #to check the dataset

animals[c(6,8)] <- c("cat", "frog")  #replacing the elements 6th and 8th 
animals

data (cars)   #loading the pre-existing dataset

data(cars)

x<- cars$speed  #assigning x axis with the values of speed
y<- cars$dist   #assigning y axis with values of dist 

#plotting the dataset
plot (x,y,las=1,
      xlab = "speed",
      ylab = "distance",
      main = "Plot between speed vs distance",
      col = "purple")


library(tidyverse)

data (starwars)

tallest_individual <- starwars[which.max(starwars$height), ]  
tallest_individual  #calling individual with the tallest height

lowest_mass_individual <- starwars[which.min(starwars$mass), ]
lowest_mass_individual   #calling individual with lowest mass




library(dplyr)
data("storms")

num_storms <- length(unique(storms$name))  #finding how many unique storms are there
num_storms

storm_counts <- table(storms$year)  #creating table for storm count
year_most_storms <- names(storm_counts[storm_counts == max(storm_counts)]) #finding year with most number of storms 
year_most_storms

highest_pressure <- max(storms$pressure, na.rm = TRUE)  #finding the max pressure of a storm in the table
storm_highest_pressure <- storms[storms$pressure == highest_pressure, c("name", "year", "pressure")]  #printing the required details like name, year and pressure
storm_highest_pressure


# Load ggplot2
library(ggplot2)

data("diamonds")
head(diamonds)

#Scatter plot of carat vs price
ggplot(diamonds, aes(x = carat, y = price))+
  geom_point()+
  labs(title = "Carat vs. Price", x = "Carat", y = "Price")

#Histogram of prices
ggplot(diamonds, aes(x = price))+
  geom_histogram(binwidth = 500, color = "brown", fill = "maroon")+
  labs(title = "Histogram of Diamond Prices", x = "Price", y = "Count")

#Boxplot of price by cut
ggplot(diamonds, aes(cut, price))+
  geom_boxplot()+
  labs(title = "Boxplot of Price by Cut", x = "Cut", y = "Price")

# Find max and min prices
max_price <- max(diamonds$price)
min_price <- min(diamonds$price)

# Add a column to classify prices as "Low price" or "High price"
diamonds$price_category <- ifelse(diamonds$price < 5000, "Low price", "High price")

# Display results
max_price
min_price
head(diamonds)

# Define a function to convert carats to grams or milligrams
carat_to_weight <- function(carat, unit = "grams") {
  if (unit == "grams") {
    return(carat * 0.2)
  } else if (unit == "milligrams") {
    return(carat * 200)
  } else {
    stop("Please choose 'grams' or 'milligrams'.")
  }
}

# Convert carat to grams for all diamonds and find the heaviest diamond
diamonds$weight_grams <- carat_to_weight(diamonds$carat, "grams")
heaviest_diamond <- max(diamonds$weight_grams)

# Display the heaviest diamond weight
heaviest_diamond






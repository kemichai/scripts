"""
...

Location: Chavannes-pres-renens, CH
Date: Feb 2023
Author: Konstantinos Michailos
"""
print("Think of a number between 0 and 100...")
# TODO: Need to change this so it does not print the guessed number in the monitor...
x = input('Enter your secret number: ')

# At the start the highest the number could be is 100 and the lowest is 0.
low = 0
high = 100
guess = int((high + low) / 2)

# Loop until we guess it correctly
while abs(guess) != int(x):
    print("--- --- --- --- --- --- --- --- --- --- --- --- ")
    print(f"Is your secret number {guess}?")
    print("--- --- --- --- --- --- --- --- --- --- --- --- ")

    print("Enter 'h' to indicate the guess is too high.")
    print("Enter 'l' to indicate the guess is too low.")
    user_inp = (input("Enter letter :"))
    if user_inp == 'h':
        high = guess
        guess = int(abs(high + low) / 2)
    elif user_inp == 'l':
        # Guess was too low. So make the current guess the lowest possible guess.
        low = guess
        guess = int(abs(high + low) / 2)

    else:
        assert False, (
                "OH NO! The program only takes as an input "
                "'h' or 'l'... ")

print("Wait I think we found your secret number... Was it:"),
print(guess)

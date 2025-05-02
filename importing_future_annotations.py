
# avoiding a heavy import at runtime 
from __future__ import annotations
from datetime import datetime 
from typing import ClassVar

class CookieBatch:
    def __init__(self, flavor):
        self.flavor = flavor

    def taste(self):
        print(f"These {self.flavor} cookies are yummy!")

# Bake one batch and try a bite 
snack = CookieBatch("chocolatechip")
snack.taste() 


class CookieBatch:
    def __init__(self, flavor: str, pieces: int = 12) -> None: 
        self.flavor = flavor 
        self.pieces = pieces 

def bake(self, minutes: int) -> None:
    print(f"Baking {self.pieces} {self.flavor} cookies for {minutes} min...") 


import os
import numpy

def classify_query_length(query):
    words_number = len(query.strip().split())
    if words_number <= 5:
        return "Short"
    elif 6 <= words_number <= 12:
        return "Medium"
    else:
        return "Long"
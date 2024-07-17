from Bio import Entrez


def get_taxid_from_accession(accession):
  # Email address required by Entrez
  Entrez.email = 'pablo.alessandro@gmail.com'
  Entrez.api_key = '86cf88e5ee6087442f57c78ed90336b99408'    
  try:
    # Search for the accession number
    search_handle = Entrez.esearch(db="nuccore", term=accession)
    search_results = Entrez.read(search_handle)
    search_handle.close()

    if not search_results['IdList']:
      return None, "No records found."
    
    print(f"Found {len(search_results['IdList'])} tax id record for accession {accession}")

    # Fetch the summary of the search results
    fetch_handle = Entrez.esummary(db="nuccore", id=search_results['IdList'][0])
    summary = Entrez.read(fetch_handle)
    fetch_handle.close()

    # Extract TaxID
    taxid = int(summary[0]['TaxId'])
    print(f"Tax id for accession {accession} is {taxid}")
    return taxid, ""

  except Exception as e:
    return None, str(e)




def main():
  # Example usage
  accession_number = "NC_045512"  # Replace with your accession number
  taxid, error = get_taxid_from_accession(accession_number)

  if error:
    print(f"Error: {error}")
  else:
    print(f"TaxID for accession {accession_number} is {taxid}")

  # with open(output_file, 'w') as file:
  #   file.write(output_content)

    
if __name__ == '__main__':
    main()

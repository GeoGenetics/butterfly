from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
import os
import time

download_dir = os.path.join("/maps/projects/caeg/scratch/dlm551/robotbutterfly/", "downloads")  

if not os.path.exists(download_dir):
    os.makedirs(download_dir)

chrome_options = webdriver.ChromeOptions()
prefs = {
    "download.default_directory": download_dir,
    "download.prompt_for_download": False,
    "download.directory_upgrade": True,
    "safebrowsing.enabled": True  
}
chrome_options.add_experimental_option("prefs", prefs)
chrome_options.add_argument("--allow-running-insecure-content")
chrome_options.add_argument("--ignore-certificate-errors")
chrome_options.add_argument("--unsafely-treat-insecure-origin-as-secure=http://dandyweb01fl.unicph.domain:5100/search") # this is the magical thing that makes it possible to download
chrome_options.add_argument("--disable-features=InsecureDownloadWarnings")
chrome_options.add_argument("--disable-web-security")  
chrome_options.add_argument("--headless") # this is needed for cron-ness 
chrome_options.add_argument("--no-sandbox")

driver = webdriver.Chrome(options=chrome_options)

driver.get("http://dandyweb01fl.unicph.domain:5100/search")

try:
    wait = WebDriverWait(driver, 10) #checkbox_smdb
    checkbox = wait.until(EC.element_to_be_clickable((By.ID, "checkbox_smdb")))
    checkbox.click()

    download_button = WebDriverWait(driver, 10).until(
        EC.element_to_be_clickable((By.XPATH, "//button[text()='Download joined table...']"))
    )
    download_button.click()
    time.sleep(10)

    # maybe this wont work in future, keep page open until download suffix isnt weird 
    download_complete = False
    while not download_complete:
        time.sleep(2) 
        download_complete = all(
            not filename.endswith(".crdownload") for filename in os.listdir(download_dir)
        )

except Exception as e:
    print("Error:", e)
finally:
    driver.quit()